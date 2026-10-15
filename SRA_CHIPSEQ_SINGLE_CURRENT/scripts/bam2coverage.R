#!/usr/bin/env Rscript
# Usage: bam2coverage.R <in.bam> <fragment_length> <out.rds> <out.bw>
# Reads are extended to <fragment_length> from their 5' end (strand-aware);
# <fragment_length> <= 0 keeps the aligned read as is.
suppressWarnings(suppressMessages({
  library(GenomicAlignments)
  library(rtracklayer)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("Usage: bam2coverage.R <in.bam> <fragment_length> <out.rds> <out.bw>")
fragment.length <- as.integer(args[2])

reads <- granges(readGAlignments(args[1]))
saveRDS(reads, file = args[3])
# extension can run past a chromosome end: trim to the chromosome lengths
if (fragment.length > 0) reads <- suppressWarnings(trim(resize(reads, fragment.length)))
export.bw(coverage(reads), args[4])
