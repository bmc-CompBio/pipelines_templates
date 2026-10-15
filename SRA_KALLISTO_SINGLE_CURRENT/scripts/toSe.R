library(tximport)
library(SummarizedExperiment)
library(biomaRt)
library(GenomicFeatures)
library(txdbmaker)

cluster.dir <- "./"


yamlDat <- yaml::yaml.load_file(paste0(cluster.dir,"config.yaml"))
annotVersion <- yamlDat$annotationVersion
organism <- yamlDat$organism

### SAMPLES TABLE
samples <- read.delim(paste0(cluster.dir,yamlDat$runtable), sep=",")
samples <- samples[!duplicated(samples$Sample.Name),]
samples <- data.frame(samples, row.names = samples$sample.id)


gtf <- paste0(cluster.dir, "genome.gtf.gz")
txdb <- makeTxDbFromGFF(gtf)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

files <- paste0(cluster.dir, samples$Sample.Name, "_out/",samples$Sample.Name,".abundance.tsv")
txi <- tximport(files, type="kallisto", tx2gene = tx2gene, ignoreTxVersion=T)


cmat <- round(txi$counts)
mode(cmat) <- "integer"
tmat <- txi$abundance
avgTxLength <- txi$length

flt <- rowSums(cmat) > 0
cmat <- cmat[flt, , drop=FALSE]
tmat <- tmat[flt, , drop=FALSE]
avgTxLength <- avgTxLength[flt, , drop=FALSE]
dim(cmat)

rownames(cmat) <- unname(sapply(rownames(cmat), function(x){strsplit(x, "\\.")[[1]][1]}, USE.NAMES = F))
rownames(tmat) <- unname(sapply(rownames(tmat), function(x){strsplit(x, "\\.")[[1]][1]}, USE.NAMES = F))
colnames(cmat) <- samples$sample.id
colnames(tmat) <- samples$sample.id

datasetOrg <- switch(organism,
  "mouse"="mmusculus",
  "human"="hsapiens",
  "drosophila"="dmelanogaster")

## get Annotation data
httr::set_config(httr::config(ssl_verifypeer = FALSE))
bm <- getBM(filters = "ensembl_gene_id",
            attributes = c("ensembl_gene_id","entrezgene_id","external_gene_name","gene_biotype","description"),
            values = rownames(cmat),
            mart =  useEnsembl(biomart = "ensembl",
                               dataset = paste0(datasetOrg,"_gene_ensembl"),
                               mirror  = "useast"))

annot <- data.frame(id=rownames(cmat))
annot$symbol <- bm$external_gene_name[match(rownames(cmat), bm$ensembl_gene_id)]
annot$entrez <- bm$entrezgene[match(rownames(cmat), bm$ensembl_gene_id)]
annot$biotype <- bm$gene_biotype[match(rownames(cmat), bm$ensembl_gene_id)]
annot$description <- bm$description[match(rownames(cmat), bm$ensembl_gene_id)]


se <- SummarizedExperiment(assays=list(counts=cmat, tpms=tmat, avgTxLength=avgTxLength), colData=samples, rowData=annot)

metadata(se)$annotationVersion <- annotVersion
metadata(se)$transcripts <- yamlDat$transcriptsFTP
metadata(se)$gtf <- yamlDat$gtfFTP
metadata(se)$strandedness <- yamlDat$strandedness
metadata(se)$organism <- yamlDat$organism 

dim(se)
saveRDS(se, file="se.rds")
