# SRA ChIP-seq (single-end)

Snakemake pipeline that takes single-end ChIP-seq runs from NCBI SRA to unique, deduplicated BAMs, MACS2 peak calls, fragment-extended coverage tracks and a MultiQC report. It is meant to be copied and adapted: for a new dataset only `config.yaml`, `SraRunTable.csv` and `peakcall_map.txt` change. It is the single-end counterpart of `SRA_CHIPSEQ_PAIRED_CURRENT`.

```
prefetch → fastq-dump → bowtie2 → unique reads only → samtools markdup -r → BAM + index
                                                          ├→ MACS2 (--nomodel --extsize <fragment_length>) → peak_calls/
                                                          ├→ reads extended to <fragment_length> → bigWig + GRanges (.rds)
                                                          └→ flagstat / stats → MultiQC
```

1. **Download** – `prefetch` (login node, needs internet) and `fastq-dump` (compute node) are separate rules. Runs with the same `Sample Name` are merged at the alignment step.
2. **Align** – `bowtie2 --no-unal`. Reads with a second-best alignment (`XS:i` tag) are discarded, so only uniquely mapping reads remain.
3. **Deduplicate** – `samtools sort | markdup -r`; the BAM is coordinate-sorted and indexed.
4. **Peaks** – `macs2 callpeak -f BAM --nomodel --extsize <fragment_length> -p <macs2_pvalue>`, with the effective genome size (non-N bases) counted from the FASTA unless given in the config.
5. **Coverage** – every read is extended to `fragment_length` bp in its 5'→3' direction (strand-aware, trimmed at chromosome ends) and written as bigWig plus a `GRanges` of the unextended reads. Coverage is raw (not normalised).
6. **QC** – `samtools flagstat/stats` and the bowtie2 logs are aggregated by MultiQC.

## Fragment length

`fragment_length` in `config.yaml` (default 200) sets how far single-end reads are extended for the coverage tracks, and is also passed to MACS2 as `--extsize`, so tracks and peak calls agree. Use the mean fragment size of the library (typically 150–300 bp for sonicated chromatin). `0` disables the extension (coverage of the sequenced read only; MACS2 then falls back to 200). Changing it reruns the coverage and MACS2 rules.

## Adapting to a new dataset

| File | What to change |
|---|---|
| `config.yaml` | `genomeFTP`: URL of the gzipped reference FASTA. `fragment_length`, `macs2_pvalue`, `macs2_broad` (`true` for broad marks). Optional: `bowtie_index` + `genome_size` (reuse an existing index), `max_spots` (test mode). |
| `SraRunTable.csv` | From the NCBI SRA Run Selector ("Metadata"). Needs the columns `Run` and `Sample Name`. |
| `peakcall_map.txt` | Tab-separated `enriched`, `control`, `name`; one row per peak call. `enriched`/`control` are `Sample Name` values from the run table. Repeating a `name` pools its samples in one MACS2 call. |

The shipped example is the human T-cell study PRJNA479884 (cFOS, JUNB, NFAT1 ChIP-seq, activated vs. resting CD4 T cells, GRCh38). In earlier versions the peak-call pairs were derived from the `chip_antibody` / `treatment` columns of that table; they are now listed explicitly in `peakcall_map.txt`, as in the paired-end pipeline.

## Installation

```bash
conda env create -f environment.yml
conda activate sra_chipseq_single
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
| `bowtie` (align + filter + dedup) | 16 | 32 GB | 480 | slim16 |
| `bam_qc` | 2 | 2 GB | 60 | slim16 |
| `bam2rangeAndcoverage` | 1 | 16 GB | 120 | slim16 |
| `macs2` | 1 | 16 GB | 120 | slim16 |
| `multiqc` | 1 | 4 GB | 30 | slim16 |

### Quick test

Set `max_spots: 10000000` in `config.yaml` (only the first 10 M spots per run are dumped), reduce `SraRunTable.csv` / `peakcall_map.txt` to one ChIP and one control sample, and point `genomeFTP` at a single chromosome FASTA. `prefetch` always downloads the complete SRA file.

## Output

| Path | Description |
|---|---|
| `peak_calls/{name}_peaks.narrowPeak` (`.broadPeak` if `macs2_broad`) | MACS2 peaks |
| `peak_calls/{name}_macs2.log` | MACS2 log |
| `{sample}_out/{sample}.bam` (+ `.bai`) | Unique, deduplicated, coordinate-sorted reads |
| `{sample}_out/{sample}_coverage.bw` | Coverage of reads extended to `fragment_length` |
| `{sample}_out/{sample}_ranges.rds` | Unextended reads as `GRanges` |
| `{sample}_out/alignstats.txt`, `markdup_stats.txt`, `.flagstat`, `.stats` | Alignment and duplicate statistics |
| `qc/multiqc_report.html` | MultiQC summary |
| `genome_size.txt`, `bowtie_index/` | Intermediate; the index is kept so it can be reused (`bowtie_index:` in the config) |

Raw FASTQs, SRA files and the genome FASTA are temporary and removed by Snakemake.

## Notes

- Duplicates are removed (`markdup -r`), not just flagged.
- License: MIT (see repository root).
