# pipelines_templates

Snakemake pipelines used directly or as templates. Each folder is self-contained and has its own README.

| Pipeline | Description | Status |
|---|---|---|
| [`SRA_CHIPSEQ_PAIRED_CURRENT`](SRA_CHIPSEQ_PAIRED_CURRENT/) | Paired-end ChIP-seq from SRA: bowtie2, dedup, MACS2, coverage, MultiQC | modernized, tested on the LMU cluster |
| [`SRA_CHIPSEQ_SINGLE_CURRENT`](SRA_CHIPSEQ_SINGLE_CURRENT/) | Single-end ChIP-seq from SRA: unique reads, MACS2, fragment-extended coverage (`fragment_length` in config), MultiQC | modernized, tested on the LMU cluster |
| [`SRA_ATAC_PAIRED_CURRENT`](SRA_ATAC_PAIRED_CURRENT/) | Paired-end ATAC-seq from SRA: fastp, bowtie2, MAPQ/dedup/mito filtering, MACS2, coverage, MultiQC | modernized, tested on the LMU cluster |
| `SRA_KALLISTO_SINGLE_CURRENT` | Single-end RNA-seq quantification (kallisto) from SRA | not yet updated |

Modernized pipelines target the LMU SLURM cluster (`profile/` directory, single conda environment per pipeline).

## Running on the LMU cluster

The pipelines are executed with Claude Code's `cluster-execution` skill: it copies a pipeline folder to the cluster, builds the conda environment, writes/uses the SLURM profile, runs Snakemake, retrieves the files declared in `rule all`, and cleans up on success. A sanitized copy of the skill is in [`docs/cluster-execution.md`](docs/cluster-execution.md).

To use it:

1. Copy it to `~/.claude/skills/cluster-execution/SKILL.md`.
2. Replace the placeholders with your own values: `<user>` (cluster account), `<cluster-login-host>` (login node, also in the `description` line), `<gpu-node>`, `<path-to>`. Adjust `/work/project/<user>` if your work directory differs.
3. In Claude Code, run `/cluster-execution <path-to-pipeline-folder>`.

Without Claude Code, the file doubles as documentation of the manual procedure (rsync, `mamba env create`, `snakemake --profile profile/`) and of cluster constraints (no internet on compute nodes, `sra-tools` pinned to 2.11.0, GPU `--gres`). The `profile/` folder of each pipeline is the SLURM profile used there.
