#!/usr/bin/env python3
"""
Compute the effective genome size (number of non-N bases) from a FASTA file.
Used to supply -g to MACS2 callpeak.

Usage:
    python get_genome_size.py genome.fa > genome_size.txt
"""

import sys


def count_effective_bases(fasta_path: str) -> int:
    valid = frozenset("ACGTacgt")
    total = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            total += sum(1 for c in line if c in valid)
    return total


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <genome.fa>")
    print(count_effective_bases(sys.argv[1]))
