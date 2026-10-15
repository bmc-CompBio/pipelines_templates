# conda env create --file sra_chipseq_paired.yaml
conda activate sra_chipseq_paired

snakemake -np
snakemake --dag | dot -Tsvg > dag.svg

## on local machine
snakemake --cores 16

## on the cluster
snakemake --profile slurm/
