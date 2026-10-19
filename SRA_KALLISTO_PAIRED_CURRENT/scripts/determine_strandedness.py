# Translate the strand bias reported by salmon (lib_format_counts.json) into the
# strandedness used for kallisto: < 0.2 reverse (ISR, dUTP), > 0.8 forward (ISF).
import json

with open(snakemake.input.json) as f:
    bias = json.load(f)["strand_mapping_bias"]

if bias < 0.2:
    strandedness = "reverse"
elif bias > 0.8:
    strandedness = "forward"
else:
    strandedness = "unstranded"

with open(snakemake.output.strandedness, "w") as f:
    f.write(strandedness + "\n")
