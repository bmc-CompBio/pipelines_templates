# SRA ChIP-seq (paired-end)

Snakemake pipeline that takes paired-end ChIP-seq runs from NCBI SRA to deduplicated BAMs, MACS2 peak calls, coverage tracks and a MultiQC report. It is meant to be copied and adapted: for a new dataset only `config.yaml`, `SraRunTable.csv` and `peakcall_map.txt` change.

```
prefetch → fastq-dump → bowtie2 → samtools (fixmate, markdup -r) → BAM + index
                                        ├→ MACS2 (BAMPE)           → peak_calls/
                                        ├→ fragments → bigWig + GRanges (.rds)
                                        └→ flagstat / stats        → MultiQC
```

1. **Download** – `prefetch` (login node, needs internet) and `fastq-dump --split-3` (compute node) are separate rules. Runs of the same `Sample Name` are merged at the alignment step.
2. **Align** – `bowtie2` (`-X 1000`, no discordant/mixed/unaligned), keeping proper pairs only.
3. **Deduplicate** – `samtools sort -n | fixmate -m | sort | markdup -r`, streamed with uncompressed BAM between stages. The final BAM is coordinate-sorted and indexed.
4. **Peaks** – `macs2 callpeak -f BAMPE`, using the effective genome size (non-N bases) counted from the FASTA (or given in the config).
5. **Coverage** – fragments (chrom, start, end) are taken from the upstream mate of each proper pair and turned into a bigWig and an R `GRanges` object with `rtracklayer`. Coverage is raw (not normalised). Fragment ends come from `TLEN`; compared with the former `bamToBed -bedpe` route this is identical for >99.99 % of fragments (the rest differ by 1–2 bp).
6. **QC** – `samtools flagstat/stats` and the bowtie2 logs are aggregated by MultiQC.

## Adapting to a new dataset

| File | What to change |
|---|---|
| `config.yaml` | `genomeFTP`: URL of the gzipped reference FASTA. `macs2_broad`: `false` for TFs / sharp marks, `true` for broad marks. Optional: `bowtie_index` + `genome_size` (reuse an existing index), `max_spots` (test mode). |
| `SraRunTable.csv` | From the NCBI SRA Run Selector ("Metadata"). Needs the columns `Run` and `Sample Name`. |
| `peakcall_map.txt` | Tab-separated `enriched`, `control`, `name`; one row per peak call. `enriched`/`control` are `Sample Name` values from the run table. |

The shipped example is the mouse p65 ChIP-seq study PRJNA388099 (10 runs, GRCm39).

## Installation

```bash
conda env create -f environment.yml
conda activate sra_chipseq_paired
```

One environment covers all rules; `--use-conda` is not used. `sra-tools` is pinned to 2.11.0 (newer versions need a newer glibc than the cluster's CentOS 7), which in turn limits `samtools` to <1.17.

## Running

```bash
snakemake -n                          # dry run
snakemake --cores 16                  # local machine (needs internet)
snakemake --profile profile/          # SLURM, see below
```

### On the BMC cluster

`profile/config.yaml` submits every rule through `sbatch` (cluster-generic executor). `prefetch` and the genome download are `localrules` and run on the login node because compute nodes have no internet. Threads and memory are set in the Snakefile, wall time in the profile.

| Rule | Threads | Memory | Time (min) | Partition |
|---|---|---|---|---|
| `prefetch`, `getGenomeFasta` | local | – | – | login node |
| `fastq_dump` | 1 | 2 GB | 180 | slim16 |
| `index_bowtie` | 16 | 32 GB | 120 | slim16 |
| `bowtie` (align + dedup) | 16 | 48 GB | 480 | slim16 |
| `bam_to_fragments` | 2 | 2 GB | 60 | slim16 |
| `bed2rangeAndcoverage` | 1 | 24 GB | 120 | slim16 |
| `bam_qc` | 2 | 2 GB | 60 | slim16 |
| `macs2` | 1 | 16 GB | 120 | slim16 |
| `multiqc` | 1 | 4 GB | 30 | slim16 |

### Quick test

Set `max_spots: 10000000` in `config.yaml` (only the first 10 M spots per run are dumped), reduce `SraRunTable.csv` / `peakcall_map.txt` to one ChIP and one input sample, and point `genomeFTP` at a single chromosome FASTA. This was used to validate the pipeline on the cluster (about 35 min wall time, dominated by `prefetch`/`fastq-dump`). `prefetch` always downloads the complete SRA file.

## Output

| Path | Description |
|---|---|
| `peak_calls/{name}_peaks.narrowPeak` (`.broadPeak` if `macs2_broad`) | MACS2 peaks |
| `peak_calls/{name}_macs2.log` | MACS2 log |
| `{sample}_out/{sample}.bam` (+ `.bai`) | Deduplicated, coordinate-sorted proper pairs |
| `{sample}_out/{sample}_coverage.bw` | Fragment coverage bigWig |
| `{sample}_out/{sample}_ranges.rds` | Fragments as `GRanges` |
| `{sample}_out/alignstats.txt`, `markdup_stats.txt`, `.flagstat`, `.stats` | Alignment and duplicate statistics |
| `qc/multiqc_report.html` | MultiQC summary |
| `genome_size.txt`, `bowtie_index/` | Intermediate; the index is kept so it can be reused (`bowtie_index:` in the config) |

Raw FASTQs, SRA files, the genome FASTA and the fragment BED are temporary and removed by Snakemake.

## Notes

- Duplicates are removed (`markdup -r`), not just flagged.
- Repeating a `name` in `peakcall_map.txt` pools its `enriched` / `control` samples in one MACS2 call.
- License: MIT (see repository root).
