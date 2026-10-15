# Input: name-sorted SAM (with header) of proper pairs.
# Keeps a pair only if both mates have MAPQ >= minq (so no orphans are left
# behind). Usage: awk -F'\t' -v minq=30 -f filter_mapq_pairs.awk
BEGIN { OFS = "\t"; have = 0 }
/^@/ { print; next }
{
    if (have && $1 == prev_name) {
        if (prev_q >= minq && $5 >= minq) { print prev_line; print $0 }
        have = 0
    } else {
        prev_line = $0; prev_name = $1; prev_q = $5; have = 1
    }
}
