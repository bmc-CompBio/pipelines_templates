# SRA ChIP-seq Pipeline (Paired-end)

Snakemake pipeline for paired-end ChIP-seq data deposited at NCBI SRA. To reuse for any other dataset, only three files need to be changed: `config.yaml`, `SraRunTable.csv`, and `peakcall_map.txt`.

## Pipeline overview

```
SRA prefetch → fastq-dump → bowtie2 → markdup → MACS2 peak calls
                                     ↘ coverage track (.bw) + GRanges (.rds)
```

1. **Download**: `prefetch` + `fastq-dump` fetch paired FASTQ files from SRA
2. **Align**: `bowtie2` aligns to the reference genome (multi-run samples are merged)
3. **Deduplicate**: `samtools fixmate` + `markdup` remove PCR duplicates
4. **Peak calling**: `macs2 callpeak -f BAMPE` using name-sorted deduplicated BAMs
5. **Coverage**: `bamToBed` + `rtracklayer` produce a bigWig track and an R GRanges object

Temporary files (`genome.fa`, raw FASTQs, intermediate BED) are deleted automatically by Snakemake after they are no longer needed. The genome size passed to MACS2 is computed from `genome.fa` by `scripts/get_genome_size.py` (counts non-N bases).

## Adapting to a new dataset

| File | What to change |
|---|---|
| `config.yaml` | `genomeFTP`: Ensembl/UCSC FTP URL for the reference genome FASTA; `macs2_broad`: `false` for TF/sharp marks, `true` for broad histone marks |
| `SraRunTable.csv` | Download from NCBI SRA Run Selector for your study; must contain columns `Run` and `Sample Name` |
| `peakcall_map.txt` | Tab-separated: `enriched`, `control`, `name` — one row per peak call, values are `Sample Name` entries from the run table |

## Installation

```bash
conda env create --file sra_chipseq_paired.yaml
conda activate sra_chipseq_paired
```

## Running

```bash
# Dry run
snakemake -np

# Visualise DAG
snakemake --dag | dot -Tsvg > dag.svg

# Local execution
snakemake --cores 16

# SLURM cluster
snakemake --profile slurm/
```

## Output files

| Path | Description |
|---|---|
| `peak_calls/{name}_peaks.narrowPeak` | MACS2 peak calls (NarrowPeak format) |
| `peak_calls/{name}_peaks.broadPeak` | MACS2 peak calls when `macs2_broad: true` |
| `peak_calls/{name}_macs2.log` | MACS2 log with FRiP and model info |
| `{sample}_out/{sample}.bam` | Name-sorted, deduplicated BAM |
| `{sample}_out/{sample}_coverage.bw` | BigWig coverage track |
| `{sample}_out/{sample}_ranges.rds` | R GRanges object (rtracklayer) |
| `{sample}_out/alignstats.txt` | bowtie2 alignment summary |
| `{sample}_out/markdup_stats.txt` | Duplicate marking statistics |
| `genome_size.txt` | Effective genome size used for MACS2 |

## Notes

- The `bowtie_index/` directory is **not** deleted automatically. Remove it manually once the pipeline has completed if disk space is needed.
- MACS2 receives the effective genome size (non-N bases) computed directly from the downloaded FASTA — no manual lookup required.
- Duplicate removal uses `samtools fixmate -m | samtools markdup -r`, operating on paired-end reads in streaming mode to avoid writing intermediate BAM files.
