# pipelines_templates

Snakemake pipelines used directly or as templates. Each folder is self-contained and has its own README.

| Pipeline | Description | Status |
|---|---|---|
| [`SRA_CHIPSEQ_PAIRED_CURRENT`](SRA_CHIPSEQ_PAIRED_CURRENT/) | Paired-end ChIP-seq from SRA: bowtie2, dedup, MACS2, coverage, MultiQC | modernized, tested on the LMU cluster |
| [`SRA_CHIPSEQ_SINGLE_CURRENT`](SRA_CHIPSEQ_SINGLE_CURRENT/) | Single-end ChIP-seq from SRA: unique reads, MACS2, fragment-extended coverage (`fragment_length` in config), MultiQC | modernized, tested on the LMU cluster |
| `SRA_KALLISTO_SINGLE_CURRENT` | Single-end RNA-seq quantification (kallisto) from SRA | not yet updated |

Modernized pipelines target the LMU SLURM cluster (`profile/` directory, single conda environment per pipeline).
