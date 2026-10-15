#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library(rtracklayer)))


args <- commandArgs(trailingOnly=TRUE)
file.id <- args[1]

bedR <- import.bed(file.id)
file.id <- strsplit(file.id,"\\.")[[1]][1]
saveRDS(bedR, file=paste(file.id,"_ranges.rds", sep=""))
cov <- coverage(bedR)
export.bw(cov, paste(file.id,"_coverage.bw", sep=""))
