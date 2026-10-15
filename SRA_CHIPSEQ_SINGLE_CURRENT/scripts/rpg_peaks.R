library(GenomicRanges)
library(rtracklayer)

BASE <- "/Volumes/work/obst/GSE116695_new"

# Load RPGs
gtf <- import(file.path(BASE, "Homo_sapiens.GRCh38.115.gtf.gz"),
              format = "gtf", feature.type = "gene",
              colnames = c("gene_name", "gene_biotype"))

ribo_gr <- gtf[grepl("^RP[LS][0-9]", gtf$gene_name) &
               gtf$gene_biotype == "protein_coding" &
               as.character(seqnames(gtf)) %in% as.character(1:22)]
ribo_gr <- ribo_gr[!duplicated(ribo_gr$gene_name)]

tss_gr <- GRanges(
  seqnames = seqnames(ribo_gr),
  ranges   = IRanges(
    start = ifelse(as.character(strand(ribo_gr)) == "+", start(ribo_gr), end(ribo_gr)),
    width = 1),
  strand = strand(ribo_gr),
  gene   = ribo_gr$gene_name)

promoter_gr <- suppressWarnings(promoters(tss_gr, upstream = 2000, downstream = 200))

# Load all peak BED files and build union
bed_files  <- list.files(file.path(BASE, "peak_calls"), pattern = "\\.bed$", full.names = TRUE)
all_peaks  <- do.call(c, lapply(bed_files, function(f) import(f, format = "BED")))
union_peaks <- reduce(all_peaks)

# Find which promoters overlap any peak
ol        <- findOverlaps(promoter_gr, union_peaks)
hit_genes <- sort(unique(tss_gr$gene[queryHits(ol)]))

cat(length(hit_genes), "RPG promoters with at least one peak (any condition):\n")
cat(paste(hit_genes, collapse = "\n"), "\n")

# Also show which comparisons each gene has a peak in
bed_names <- sub("_peaks.*", "", sub("\\.bed$", "", basename(bed_files)))
detail <- lapply(seq_along(bed_files), function(i) {
  pk <- import(bed_files[i], format = "BED")
  ol2 <- findOverlaps(promoter_gr, pk)
  genes <- tss_gr$gene[queryHits(ol2)]
  data.frame(gene = genes, comparison = bed_names[i], stringsAsFactors = FALSE)
})
detail_df <- do.call(rbind, detail)
detail_df <- detail_df[detail_df$gene %in% hit_genes, ]

agg <- aggregate(comparison ~ gene, data = detail_df,
                 FUN = function(x) paste(sort(unique(x)), collapse = "; "))
agg <- agg[order(agg$gene), ]

cat("\nDetailed breakdown:\n")
print(agg, row.names = FALSE)
