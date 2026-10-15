#!/usr/bin/env Rscript
# Usage: bed2rAc.R <fragments.bed> <ranges.rds> <coverage.bw>
suppressWarnings(suppressMessages(library(rtracklayer)))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: bed2rAc.R <in.bed> <out.rds> <out.bw>")

bedR <- import.bed(args[1])
saveRDS(bedR, file = args[2])
export.bw(coverage(bedR), args[3])
