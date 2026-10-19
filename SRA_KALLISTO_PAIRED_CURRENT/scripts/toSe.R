# Build a SummarizedExperiment (counts, TPM, avg. transcript length) from the
# kallisto abundances with tximport; gene annotation from the Ensembl GTF and FTP mappings.
suppressPackageStartupMessages({
  library(tximport)
  library(SummarizedExperiment)
  library(GenomicFeatures)
})

cfg <- yaml::yaml.load_file("config.yaml")

## samples: one row per Sample Name
samples <- read.delim(cfg$runtable, sep = ",", check.names = FALSE)
samples <- samples[!duplicated(samples[["Sample Name"]]), ]
sample.names <- samples[["Sample Name"]]
samples <- data.frame(samples, row.names = sample.names, check.names = FALSE)

## transcript -> gene
makeTx <- if (requireNamespace("txdbmaker", quietly = TRUE)) {
  txdbmaker::makeTxDbFromGFF
} else {
  GenomicFeatures::makeTxDbFromGFF
}
txdb <- makeTx("genome.gtf.gz")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

files <- paste0(sample.names, "_out/", sample.names, ".abundance.tsv")
names(files) <- sample.names
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)

cmat <- round(txi$counts)
mode(cmat) <- "integer"
tmat <- txi$abundance
avgTxLength <- txi$length

flt <- rowSums(cmat) > 0
cmat <- cmat[flt, , drop = FALSE]
tmat <- tmat[flt, , drop = FALSE]
avgTxLength <- avgTxLength[flt, , drop = FALSE]
dim(cmat)

strip.version <- function(x) sub("\\..*$", "", x)
rownames(cmat) <- strip.version(rownames(cmat))
rownames(tmat) <- strip.version(rownames(tmat))
rownames(avgTxLength) <- strip.version(rownames(avgTxLength))

## gene annotation, all from the Ensembl release given in config.yaml (no network access needed):
# symbol and biotype from the GTF, entrez id from the Ensembl FTP "entrez" mapping
gtf.genes <- rtracklayer::import("genome.gtf.gz", feature.type = "gene")
g <- mcols(gtf.genes)
gtf.annot <- data.frame(id = g$gene_id,
                        symbol = if (is.null(g$gene_name)) NA else g$gene_name,
                        biotype = if (is.null(g$gene_biotype)) NA else g$gene_biotype)
m <- match(rownames(cmat), gtf.annot$id)

# EntrezGene rows only ("EntrezGene_trans_name" are name-based guesses); a gene with
# several entrez ids gets the one that is listed for most of its transcripts
ez <- read.delim("entrez.tsv.gz", colClasses = "character")
ez <- ez[ez$db_name == "EntrezGene", c("gene_stable_id", "xref")]
if (nrow(ez) > 0) {
  ez <- as.data.frame(table(ez$gene_stable_id, ez$xref), stringsAsFactors = FALSE)
  ez <- ez[ez$Freq > 0, ]
  ez <- ez[order(-ez$Freq, as.numeric(ez$Var2)), ]
  ez <- ez[!duplicated(ez$Var1), ]
}
me <- match(rownames(cmat), ez$Var1)

annot <- data.frame(id = rownames(cmat),
                    symbol = gtf.annot$symbol[m],
                    entrez = as.integer(ez$Var2[me]),
                    biotype = gtf.annot$biotype[m])
rownames(annot) <- annot$id
message(sum(!is.na(annot$entrez)), " of ", nrow(annot), " genes have an entrez id")

se <- SummarizedExperiment(assays = list(counts = cmat, tpms = tmat, avgTxLength = avgTxLength),
                           colData = samples, rowData = annot)

metadata(se)$annotationVersion <- cfg$annotationVersion
metadata(se)$transcripts <- cfg$transcriptsFTP
metadata(se)$gtf <- cfg$gtfFTP
metadata(se)$strandedness <- if (file.exists("results/strandedness.txt")) readLines("results/strandedness.txt") else cfg$strandedness
metadata(se)$organism <- cfg$organism

dim(se)
saveRDS(se, file = "se.rds")
