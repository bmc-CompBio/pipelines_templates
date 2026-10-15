# SRA ATAC-seq (paired-end)

Snakemake pipeline that takes paired-end ATAC-seq runs from NCBI SRA to filtered, deduplicated BAMs without mitochondrial reads, MACS2 peak calls, fragment coverage tracks and a MultiQC report. It is meant to be copied and adapted: for a new dataset only `config.yaml` and `SraRunTable.csv` change.

```
prefetch → fastq-dump → fastp → bowtie2 → proper pairs, MAPQ filter → fixmate → markdup -r
                                                                          ↓
                                              remove mitochondrial reads → BAM + index
                                                  ├→ MACS2 (BAMPE, no control) → peak_calls/
                                                  ├→ fragments → bigWig + GRanges (.rds)
                                                  └→ flagstat / stats / idxstats → MultiQC
```

1. **Download** – `prefetch` (login node, needs internet) and `fastq-dump --split-3` (compute node) are separate rules. Runs of the same `Sample Name` are merged at the alignment step.
2. **Trim** – `fastp` removes adapters (auto-detected from the read pairs, which catches Nextera read-through) and very short reads.
3. **Align** – `bowtie2 -X <bowtie_maxins>` (default 2000, so that multi-nucleosome fragments are kept), no discordant/mixed/unaligned pairs.
4. **Filter** – only proper pairs where **both** mates have MAPQ ≥ `min_mapq` (default 30) are kept, so no orphan mates are left behind. Duplicates are removed with `samtools fixmate -m | markdup -r` (replaces picard; no 128 GB heap needed and no second conda environment).
5. **Mitochondrial reads** – reads on `mito_name` are dropped from the final BAM. Their share is visible in the MultiQC `idxstats` plot.
6. **Peaks** – `macs2 callpeak -f BAMPE --nomodel --keep-dup all -p <macs2_pvalue>` on the fragments (duplicates are already removed), per sample without control. Genome size is counted from the FASTA unless given.
7. **Coverage** – fragments (chrom, start, end) are taken from the upstream mate of each proper pair and turned into a bigWig and a `GRanges` object with `rtracklayer`. Coverage is raw (not normalised).
8. **QC** – fastp, bowtie2, `samtools flagstat/stats` (including the insert-size histogram) and `idxstats` are aggregated by MultiQC.

## Adapting to a new dataset

| File | What to change |
|---|---|
| `config.yaml` | `genomeFTP`: URL of the gzipped reference FASTA. `mito_name`: name of the mitochondrial sequence in that FASTA (`MT` for most Ensembl genomes, `mitochondrion_genome` for *Drosophila*; empty keeps everything). Optional: `bowtie_maxins`, `min_mapq`, `macs2_pvalue`, `macs2_broad`, `bowtie_index` + `genome_size` (reuse an existing index), `max_spots` (test mode). |
| `SraRunTable.csv` | From the NCBI SRA Run Selector ("Metadata"). Needs the columns `Run` and `Sample Name`. |

The shipped example is the *Drosophila melanogaster* S2-cell ATAC-seq study PRJNA400453 / GSE103177 (2 samples, BDGP6.54).

Not included (add if needed): Tn5 cut-site shifting (+4/−5), blacklist filtering, TSS enrichment.

## Installation

```bash
conda env create -f environment.yml
conda activate sra_atac_paired
```

One environment covers all rules; `--use-conda` is not used. `sra-tools` is pinned to 2.11.0 (newer versions need a newer glibc than the cluster's CentOS 7), which in turn limits `samtools` to <1.17.

## Running

```bash
snakemake -n                          # dry run
snakemake --cores 16                  # local machine (needs internet)
snakemake --profile profile/          # SLURM, see below
```

### On the LMU cluster

`profile/config.yaml` submits every rule through `sbatch` (cluster-generic executor). `prefetch` and the genome download are `localrules` and run on the login node because compute nodes have no internet. Threads and memory are set in the Snakefile, wall time in the profile.

| Rule | Threads | Memory | Time (min) | Partition |
|---|---|---|---|---|
| `prefetch`, `getGenomeFasta` | local | – | – | login node |
| `fastq_dump` | 1 | 2 GB | 180 | slim16 |
| `fastp` | 8 | 4 GB | 120 | slim16 |
| `index_bowtie` | 16 | 32 GB | 120 | slim16 |
| `bowtie` (align + filter + dedup) | 16 | 48 GB | 480 | slim16 |
| `remove_mito` | 4 | 4 GB | 60 | slim16 |
| `bam_to_fragments` | 2 | 2 GB | 60 | slim16 |
| `bed2rangeAndcoverage` | 1 | 24 GB | 120 | slim16 |
| `bam_qc` | 2 | 2 GB | 60 | slim16 |
| `macs2` | 1 | 16 GB | 120 | slim16 |
| `multiqc` | 1 | 4 GB | 30 | slim16 |

### Quick test

Set `max_spots: 10000000` in `config.yaml` (only the first 10 M spots per run are dumped). `prefetch` always downloads the complete SRA file.

## Output

| Path | Description |
|---|---|
| `peak_calls/{sample}_peaks.narrowPeak` (`.broadPeak` if `macs2_broad`) | MACS2 peaks |
| `peak_calls/{sample}_macs2.log` | MACS2 log |
| `{sample}_out/{sample}.bam` (+ `.bai`) | Deduplicated proper pairs without mitochondrial reads, coordinate-sorted |
| `{sample}_out/{sample}_coverage.bw` | Fragment coverage bigWig |
| `{sample}_out/{sample}_ranges.rds` | Fragments as `GRanges` |
| `{sample}_out/alignstats.txt`, `markdup_stats.txt`, `.flagstat`, `.stats`, `.idxstats` | Alignment, duplicate and per-chromosome statistics (`.idxstats` is from the BAM *before* mitochondrial removal) |
| `qc/multiqc_report.html`, `qc/fastp/` | MultiQC summary and fastp reports |
| `genome_size.txt`, `bowtie_index/` | Intermediate; the index is kept so it can be reused (`bowtie_index:` in the config) |

Raw FASTQs, trimmed FASTQs, SRA files, the genome FASTA, the unfiltered BAM and the fragment BED are temporary and removed by Snakemake.

## Notes

- Earlier versions of this pipeline called MACS2 with `-g mm` regardless of the organism; the genome size is now derived from the reference FASTA.
- Sample names are the `Sample Name` column (GSM accessions); earlier versions used `Experiment` (SRX).
- License: MIT (see repository root).
