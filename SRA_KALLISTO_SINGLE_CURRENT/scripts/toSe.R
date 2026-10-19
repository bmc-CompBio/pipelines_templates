# Build a SummarizedExperiment (counts, TPM, avg. transcript length) from the
# kallisto abundances with tximport; gene annotation from biomaRt.
suppressPackageStartupMessages({
  library(tximport)
  library(SummarizedExperiment)
  library(biomaRt)
  library(GenomicFeatures)
})

cfg <- yaml::yaml.load_file("config.yaml")

## samples: one row per Sample Name
samples <- read.delim(cfg$runtable, sep = ",", check.names = FALSE)
samples <- samples[!duplicated(samples[["Sample Name"]]), ]
sample.names <- samples[["Sample Name"]]
# columns of the SummarizedExperiment: "sample.id" if the run table has it
# (e.g. a descriptive name), else the Sample Name
col.names <- if ("sample.id" %in% colnames(samples)) samples[["sample.id"]] else sample.names
samples <- data.frame(samples, row.names = col.names, check.names = FALSE)

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
names(files) <- col.names
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

## gene annotation
# symbol and biotype come from the GTF (exact for the chosen release, no network);
# entrez id and description are added from biomaRt when it is reachable
gtf.genes <- rtracklayer::import("genome.gtf.gz", feature.type = "gene")
g <- mcols(gtf.genes)
gtf.annot <- data.frame(id = g$gene_id,
                        symbol = if (is.null(g$gene_name)) NA else g$gene_name,
                        biotype = if (is.null(g$gene_biotype)) NA else g$gene_biotype)
m <- match(rownames(cmat), gtf.annot$id)
annot <- data.frame(id = rownames(cmat),
                    symbol = gtf.annot$symbol[m],
                    entrez = NA_integer_,
                    biotype = gtf.annot$biotype[m],
                    description = NA_character_)
rownames(annot) <- annot$id

datasets <- c(human = "hsapiens_gene_ensembl", mouse = "mmusculus_gene_ensembl",
              rat = "rnorvegicus_gene_ensembl", drosophila = "dmelanogaster_gene_ensembl",
              zebrafish = "drerio_gene_ensembl")
dataset <- if (!is.null(cfg$biomart_dataset)) cfg$biomart_dataset else datasets[cfg$organism]
use.biomart <- is.null(cfg$biomart) || isTRUE(cfg$biomart)
if (use.biomart && (is.null(dataset) || is.na(dataset))) {
  warning("Unknown organism '", cfg$organism, "': set biomart_dataset in config.yaml; ",
          "entrez/description left empty")
  use.biomart <- FALSE
}

if (use.biomart) {
  httr::set_config(httr::config(ssl_verifypeer = FALSE))
  host <- if (is.null(cfg$biomart_host)) "https://useast.ensembl.org" else cfg$biomart_host
  bm <- NULL
  for (attempt in 1:3) {
    bm <- tryCatch({
      mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = unname(dataset), host = host)
      getBM(filters = "ensembl_gene_id",
            attributes = c("ensembl_gene_id", "entrezgene_id", "description"),
            values = rownames(cmat), mart = mart)
    }, error = function(e) {
      message("biomaRt attempt ", attempt, " failed: ", conditionMessage(e))
      NULL
    })
    if (!is.null(bm)) break
    Sys.sleep(10)
  }
  if (is.null(bm)) {
    warning("biomaRt not reachable; entrez and description are NA")
  } else {
    bm <- bm[!duplicated(bm$ensembl_gene_id), ]
    mb <- match(rownames(cmat), bm$ensembl_gene_id)
    annot$entrez <- bm$entrezgene_id[mb]
    annot$description <- bm$description[mb]
  }
}

se <- SummarizedExperiment(assays = list(counts = cmat, tpms = tmat, avgTxLength = avgTxLength),
                           colData = samples, rowData = annot)

metadata(se)$annotationVersion <- cfg$annotationVersion
metadata(se)$transcripts <- cfg$transcriptsFTP
metadata(se)$gtf <- cfg$gtfFTP
metadata(se)$strandedness <- if (file.exists("results/strandedness.txt")) readLines("results/strandedness.txt") else cfg$strandedness
metadata(se)$organism <- cfg$organism

dim(se)
saveRDS(se, file = "se.rds")
