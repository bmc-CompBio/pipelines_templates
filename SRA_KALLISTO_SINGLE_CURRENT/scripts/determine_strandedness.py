# Translate the library type detected by salmon (--libType A) into the
# strandedness used for kallisto. Supports salmon >= 1.10 (meta_info.json,
# "library_types": SR / SF / U) and the strand bias of lib_format_counts.json.
import json

with open(snakemake.input.json) as f:
    data = json.load(f)

if "library_types" in data:
    types = data["library_types"]
    if any(t in ("ISR", "SR") for t in types):
        strandedness = "reverse"
    elif any(t in ("ISF", "SF") for t in types):
        strandedness = "forward"
    else:
        strandedness = "unstranded"
elif "strand_mapping_bias" in data:
    bias = data["strand_mapping_bias"]
    strandedness = "reverse" if bias < 0.2 else "forward" if bias > 0.8 else "unstranded"
else:
    strandedness = "unstranded"

with open(snakemake.output.strandedness, "w") as f:
    f.write(strandedness + "\n")
