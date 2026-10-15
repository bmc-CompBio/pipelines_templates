# Create/activate conda environment (one-time setup):
# mamba env create --name ChIPSeq_SRA_single --file sra_chipseq_single.yaml
# conda activate ChIPSeq_SRA_single

# ── Full dry-run ───────────────────────────────────────────────────────────────
snakemake -np

# ── DAG inspection ─────────────────────────────────────────────────────────────
snakemake --dag | dot -Tsvg > dag.svg

# ── One-sample local test: cFOS activated (anti-CD3 + anti-CD28) vs no-activation ──
# Targets: peak_calls/cFOS__GSM4106151_vs_GSM4106150.bed
#          GSM4106151_out/GSM4106151_coverage.bw  (experimental)
#          GSM4106150_out/GSM4106150_coverage.bw  (control)
snakemake -np peak_calls/cFOS__GSM4106151_vs_GSM4106150.bed \
              GSM4106151_out/GSM4106151_coverage.bw \
              GSM4106150_out/GSM4106150_coverage.bw

# Run for real (local, all cores):
# snakemake -j 8 peak_calls/cFOS__GSM4106151_vs_GSM4106150.bed \
#               GSM4106151_out/GSM4106151_coverage.bw \
#               GSM4106150_out/GSM4106150_coverage.bw

# ── Full local run ─────────────────────────────────────────────────────────────
# snakemake -j 8
