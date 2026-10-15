library(ShortRead)
library(rtracklayer)

# 1 - file.id(s)
args = commandArgs(trailingOnly=TRUE)
file.name <- args[1]

fragment.length <- 200

gr <- readGAlignments(file.name)
grs <- as(gr, "GRanges")
grsr <- resize(grs, fragment.length)
saveRDS(grs, file=sub(".bam","_ranges.rds", file.name))

covs <- coverage(grsr)

saveRDS(covs, file=sub(".bam","_coverage.rds", file.name))
export.bw(covs, sub(".bam", "_coverage.bw", file.name))
