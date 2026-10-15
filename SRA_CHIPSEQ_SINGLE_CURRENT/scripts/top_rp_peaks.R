library(GenomicRanges)
library(rtracklayer)

BASE <- "/Volumes/work/obst/GSE116695_new"

# Load RPGs — strict structural RP genes only (exclude kinases / binding proteins)
gtf <- import(file.path(BASE, "Homo_sapiens.GRCh38.115.gtf.gz"),
              format = "gtf", feature.type = "gene",
              colnames = c("gene_name", "gene_biotype"))

ribo_gr <- gtf[grepl("^RP[LS][0-9]", gtf$gene_name) &
               gtf$gene_biotype == "protein_coding" &
               as.character(seqnames(gtf)) %in% as.character(1:22)]
ribo_gr <- ribo_gr[!duplicated(ribo_gr$gene_name)]

# Exclude non-structural entries
exclude <- c("RPS6KA2","RPS6KA5","RPS6KB1","RPS6KB2","RPS6KC1","RPS19BP1")
ribo_gr <- ribo_gr[!ribo_gr$gene_name %in% exclude]

tss_gr <- GRanges(
  seqnames = seqnames(ribo_gr),
  ranges   = IRanges(
    start = ifelse(as.character(strand(ribo_gr)) == "+", start(ribo_gr), end(ribo_gr)),
    width = 1),
  strand = strand(ribo_gr),
  gene   = ribo_gr$gene_name)

promoter_gr <- suppressWarnings(promoters(tss_gr, upstream = 2000, downstream = 200))

# Load narrowPeak files (have signalValue, pValue, qValue)
np_files <- list.files(file.path(BASE, "peak_calls"),
                       pattern = "_peaks\\.narrowPeak$", full.names = TRUE)

peaks_list <- lapply(np_files, function(f) {
  pk <- import(f, format = "narrowPeak")
  # derive comparison label from filename
  label <- sub("_peaks\\.narrowPeak$", "", basename(f))
  pk$comparison <- label
  pk
})

all_peaks <- do.call(c, peaks_list)

# Find peaks overlapping any RPG promoter
ol <- findOverlaps(all_peaks, promoter_gr)
hit_peaks  <- all_peaks[queryHits(ol)]
hit_genes  <- tss_gr$gene[subjectHits(ol)]
hit_peaks$gene <- hit_genes

# Build result data frame, rank by signalValue
df <- data.frame(
  gene        = hit_peaks$gene,
  comparison  = hit_peaks$comparison,
  chr         = as.character(seqnames(hit_peaks)),
  start       = start(hit_peaks),
  end         = end(hit_peaks),
  signalValue = hit_peaks$signalValue,
  pValue      = hit_peaks$pValue,
  qValue      = hit_peaks$qValue,
  stringsAsFactors = FALSE
)

# Remove duplicate peak-gene pairs (same peak overlapping same promoter)
df <- df[!duplicated(paste(df$comparison, df$chr, df$start)), ]
df <- df[order(-df$signalValue), ]

cat("Top 10 peaks in RPL/RPS promoters (ranked by signal value):\n\n")
print(head(df, 10), row.names = FALSE)
