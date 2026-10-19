# SRA RNA-seq quantification (single-end, kallisto)

Snakemake pipeline that takes single-end RNA-seq runs from NCBI SRA to transcript abundances (kallisto) and a gene-level `SummarizedExperiment` (`se.rds`, via tximport and Ensembl annotation) plus a MultiQC report. It is the single-end counterpart of `SRA_KALLISTO_PAIRED_CURRENT`. It is meant to be copied and adapted: for a new dataset only `config.yaml` and `SraRunTable.csv` change.

```
prefetch → fastq-dump → kallisto quant (-b 30) → tximport + Ensembl annotation → se.rds
                 │             ↑                                   
                 │      strandedness.txt ← salmon on 100k reads (if strandedness: auto)
                 └→ kallisto logs → MultiQC
```

1. **Download** – `prefetch` (login node, needs internet) and `fastq-dump` (compute node) are separate rules. Runs of the same `Sample Name` (e.g. several lanes) are quantified together in one kallisto call.
2. **Annotation** – the Ensembl GTF and cDNA FASTA (plus, optionally, the ncRNA FASTA, `ncrnaFTP`) from `config.yaml` are downloaded on the login node; the kallisto index is built once and removed afterwards.
3. **Strandedness** – `strandedness: auto` runs salmon (`--libType A`) on the first 100k reads of the first run and converts the detected library type to `forward` / `reverse` / `unstranded` for kallisto (`--fr-stranded` / `--rf-stranded` / none; for single-end reads salmon reports SF / SR / U). Set `strandedness` to `none`, `forward` or `reverse` to skip salmon and its index.
4. **Quantify** – `kallisto quant --single -l <frag_mean> -s <frag_sd>` with 30 bootstraps (`bootstraps`), per sample. Single-end reads carry no fragment length information, so `frag_mean` and `frag_sd` (bp, from the library prep) must be given; they also affect the effective transcript lengths and therefore TPM.
5. **Summarise** – `scripts/toSe.R` imports the abundances with tximport (transcripts → genes via the GTF), drops genes without counts and adds annotation, all from the Ensembl release in `config.yaml` (no BioMart, no network access in this step): gene symbol and biotype from the GTF, entrez id from the Ensembl FTP file `tsv/<species>/<assembly>.<release>.entrez.tsv.gz` (URL derived from `gtfFTP`, or set `entrezFTP`). If a gene maps to several entrez ids, the one shared by most of its transcripts is used; genes without an Ensembl–Entrez link (many lncRNAs, pseudogenes, transposons) get `NA`.
6. **QC** – MultiQC summarises the kallisto logs (pseudoalignment rate, fragment length).

## Adapting to a new dataset

| File | What to change |
|---|---|
| `config.yaml` | `gtfFTP`, `transcriptsFTP`, `annotationVersion`: Ensembl GTF, cDNA FASTA and release number. `organism`: free-text label stored in the metadata of `se.rds`. `strandedness`, `frag_mean`, `frag_sd`, `bootstraps`, `ncrnaFTP`. Optional: `entrezFTP`, `max_spots` (test mode). |
| `SraRunTable.csv` | From the NCBI SRA Run Selector ("Metadata"). Needs the columns `Run` and `Sample Name`, `LibraryLayout` must be `SINGLE`. All columns are kept as `colData` of `se.rds`; if a column `sample.id` exists it is used to name the columns of `se.rds`, otherwise `Sample Name`. |

The shipped example is the *Drosophila melanogaster* DGRP RNA-seq study PRJNA483441 (957 runs, BDGP6.54, Ensembl 115; some samples have two runs). For a trial run reduce the table to a few samples.

## Installation

```bash
conda env create -f environment.yml
conda activate sra_kallisto_single
```

One environment covers all rules; `--use-conda` is not used. Two version limits come from `sra-tools`, which is pinned to 2.11.0 (newer versions need a newer glibc than the cluster's CentOS 7): `kallisto` is therefore pinned to 0.48.0 (0.50 and later need a newer HDF5), and R/Bioconductor stay at R 4.2 / Bioconductor 3.16 (no `txdbmaker`; `GenomicFeatures` provides `makeTxDbFromGFF`, and `toSe.R` uses `txdbmaker` automatically if it is installed).

## Running

```bash
snakemake -n                          # dry run
snakemake --cores 16                  # local machine (needs internet)
snakemake --profile profile/          # SLURM, see below
```

### On the BMC cluster

`profile/config.yaml` submits every rule through `sbatch` (cluster-generic executor). `prefetch`, the annotation download, `create_se` and `set_strandedness` are `localrules` and run on the login node. Threads and memory are set in the Snakefile, wall time in the profile.

| Rule | Threads | Memory | Time (min) | Partition |
|---|---|---|---|---|
| `prefetch`, `get_index_files`, `create_se` | local | – | – | login node |
| `fastq_dump` | 1 | 2 GB | 180 | slim16 |
| `kallisto_index` | 1 | 8 GB | 60 | slim16 |
| `subset_reads` | 1 | 1 GB | 30 | slim16 |
| `salmon_index` | 16 | 32 GB | 60 | slim16 |
| `salmon_alignment` | 8 | 16 GB | 60 | slim16 |
| `kallisto` | 16 | 8 GB | 240 | slim16 |
| `multiqc` | 1 | 4 GB | 30 | slim16 |

### Quick test

Set `max_spots: 2000000` in `config.yaml` (only the first 2 M spots per run are dumped) and reduce `SraRunTable.csv` to two samples (the Sample Name of a multi-run sample is a good check). `prefetch` always downloads the complete SRA file.

## Output

| Path | Description |
|---|---|
| `se.rds` | `SummarizedExperiment`: assays `counts`, `tpms`, `avgTxLength`; `rowData` gene annotation (id, symbol, entrez, biotype); `colData` SRA metadata; `metadata` with annotation sources and strandedness |
| `{sample}_out/{sample}.abundance.tsv`, `.abundance.h5` | kallisto transcript abundances (with bootstraps in the `.h5`) |
| `{sample}_out/{sample}.kallisto.log`, `run_info.json` | kallisto log and run information |
| `results/strandedness.txt` | Strandedness used for kallisto |
| `qc/multiqc_report.html` | MultiQC summary (kallisto pseudoalignment rate and fragment length per run) |
| `salmon_quant/` | Strandedness detection output (only with `strandedness: auto`) |

Raw FASTQs, SRA files, the kallisto and salmon indices and the read subsets are temporary and removed by Snakemake. `genome.gtf.gz` and `transcripts.fa.gz` are kept.

## Notes

- Raw FASTQs are not kept; the SRA files are removed after `fastq-dump`.
- Earlier versions always ran salmon (the `strandedness` option was ignored) and named the columns of `se.rds` through a hard-wired `sample.id` column; the latter is now optional.
- Earlier versions queried Ensembl BioMart for entrez id and description. BioMart is being discontinued and was unreliable, so it was replaced by the Ensembl FTP files; the `description` column was dropped.
- License: MIT (see repository root).
