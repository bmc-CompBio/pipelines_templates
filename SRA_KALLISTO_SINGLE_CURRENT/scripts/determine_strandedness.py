# scripts/determine_strandedness.py
# Supports salmon >= 1.10 (reads meta_info.json library_types)
# and older salmon (reads lib_format_counts.json strand_mapping_bias)

import json

def main():
    with open(snakemake.input.json) as f:
        data = json.load(f)

    # salmon >= 1.10: meta_info.json has library_types list e.g. ["ISR"], ["ISF"], ["U"]
    if "library_types" in data:
        lib_types = data["library_types"]
        # ISR = inward, stranded, reverse  -> reverse
        # ISF = inward, stranded, forward  -> forward
        # U / IU = unstranded
        if any(t in ("ISR", "SR") for t in lib_types):
            strandedness = "reverse"
        elif any(t in ("ISF", "SF") for t in lib_types):
            strandedness = "forward"
        else:
            strandedness = "unstranded"
    # older salmon: lib_format_counts.json has strand_mapping_bias
    elif "strand_mapping_bias" in data:
        bias = data["strand_mapping_bias"]
        if bias < 0.2:
            strandedness = "reverse"
        elif bias > 0.8:
            strandedness = "forward"
        else:
            strandedness = "unstranded"
    else:
        strandedness = "unstranded"

    with open(snakemake.output.strandedness, "w") as f:
        f.write(strandedness)

if __name__ == "__main__":
    main()
