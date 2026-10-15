suppressMessages({
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
  if (!requireNamespace("GenomicFeatures", quietly=TRUE))
    BiocManager::install("GenomicFeatures", ask=FALSE, update=FALSE)
  for (pkg in c("ggplot2", "dplyr", "patchwork")) {
    if (!requireNamespace(pkg, quietly=TRUE))
      install.packages(pkg, repos="https://cloud.r-project.org", quiet=TRUE)
  }
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

setwd("/Volumes/work/obst/GSE116695_new")

# ── 1. Ribosomal protein gene TSSs from Ensembl GTF (release 115) ────────────
gtf_gz <- "Homo_sapiens.GRCh38.115.gtf.gz"
if (!file.exists(gtf_gz)) {
  message("Downloading Ensembl GRCh38 release-115 GTF...")
  download.file(
    "https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz",
    destfile = gtf_gz, method = "wget", quiet = FALSE
  )
}

message("Parsing GTF for ribosomal protein genes...")
gtf <- import(gtf_gz, format = "gtf", feature.type = "gene",
              colnames = c("gene_name", "gene_biotype"))

# Filter: protein-coding RP[LS]<digit> on standard autosomes
ribo_gr <- gtf[grepl("^RP[LS][0-9]", gtf$gene_name) &
               gtf$gene_biotype == "protein_coding" &
               as.character(seqnames(gtf)) %in% as.character(1:22)]
ribo_gr <- ribo_gr[!duplicated(ribo_gr$gene_name)]
message(sprintf("Found %d ribosomal protein genes", length(ribo_gr)))

# TSS: 5' end of each gene
tss_gr <- GRanges(
  seqnames = seqnames(ribo_gr),
  ranges   = IRanges(
    start = ifelse(as.character(strand(ribo_gr)) == "+",
                   start(ribo_gr), end(ribo_gr)),
    width = 1
  ),
  strand = strand(ribo_gr),
  gene   = ribo_gr$gene_name
)

# ── 2. Promoter windows: TSS -2000 / +200 ────────────────────────────────────
promoter_gr <- suppressWarnings(promoters(tss_gr, upstream = 2000, downstream = 200))

# ── 3. Mean bigWig coverage over each RPG promoter (all 3 conditions) ────────
bw_files <- c(
  "No activation" = "GSM4106150_out/GSM4106150_coverage.bw",
  "Act. CD3+CD28"  = "GSM4106151_out/GSM4106151_coverage.bw",
  "Act. CD3"       = "GSM4106152_out/GSM4106152_coverage.bw"
)

mean_cov_in_region <- function(bw_file, regions) {
  sapply(seq_along(regions), function(i) {
    reg <- regions[i]
    gr  <- tryCatch(import(bw_file, which = reg, format = "BigWig"),
                    error = function(e) GRanges())
    if (length(gr) == 0) return(0)
    sum(as.numeric(width(gr)) * gr$score) / width(reg)
  })
}

message("Computing mean coverage over RPG promoters...")
cov_mat <- sapply(bw_files, mean_cov_in_region, regions = promoter_gr)
rownames(cov_mat) <- tss_gr$gene

# Fold change: max activated / (no-activation + pseudocount)
fc <- apply(cov_mat[, c("Act. CD3+CD28", "Act. CD3"), drop = FALSE], 1, max) /
      (cov_mat[, "No activation"] + 0.01)

act_signal <- apply(cov_mat[, c("Act. CD3+CD28", "Act. CD3"), drop = FALSE], 1, max)

cat("\nTop 10 RPGs by activated cFOS signal:\n")
top10 <- head(sort(act_signal, decreasing = TRUE), 10)
print(round(data.frame(
  act_max    = top10,
  fc_vs_noact = fc[names(top10)],
  no_act     = cov_mat[names(top10), "No activation"]
), 3))

# ── 4. Select 5 to plot: highest activated signal ────────────────────────────
sel_genes <- names(head(sort(act_signal, decreasing = TRUE), 5))
sel_idx   <- match(sel_genes, tss_gr$gene)
message(sprintf("\nPlotting: %s", paste(sel_genes, collapse = ", ")))

# ── 5. Per-base coverage extraction for a window ─────────────────────────────
get_cov <- function(bw_file, region, label) {
  raw  <- tryCatch(import(bw_file, which = region, format = "BigWig"),
                   error = function(e) GRanges())
  full <- data.frame(pos = seq(start(region), end(region)))
  if (length(raw) == 0) {
    full$score  <- 0
    full$sample <- label
    return(full)
  }
  # Expand each constant-score range to per-base rows
  df <- do.call(rbind, lapply(seq_along(raw), function(i)
    data.frame(pos = seq(start(raw[i]), end(raw[i])), score = raw$score[i])
  ))
  full <- left_join(full, df, by = "pos") %>%
    mutate(score = ifelse(is.na(score), 0, score), sample = label)
  full
}

# ── 6. One panel per gene ─────────────────────────────────────────────────────
plot_promoter <- function(idx) {
  prom  <- promoter_gr[idx]
  gene  <- tss_gr$gene[idx]
  tss_x <- ifelse(as.character(strand(tss_gr[idx])) == "+",
                  start(tss_gr[idx]), end(tss_gr[idx]))

  cov <- bind_rows(
    get_cov(bw_files["No activation"], prom, "No activation"),
    get_cov(bw_files["Act. CD3+CD28"],  prom, "Act. CD3+CD28"),
    get_cov(bw_files["Act. CD3"],       prom, "Act. CD3")
  ) %>% mutate(sample = factor(sample,
                  levels = c("No activation", "Act. CD3", "Act. CD3+CD28")))

  ggplot(cov, aes(pos, score, fill = sample)) +
    geom_area(alpha = 0.85, show.legend = FALSE) +
    facet_wrap(~sample, ncol = 1, scales = "free_y", strip.position = "right") +
    scale_fill_manual(values = c(
      "No activation" = "#999999",
      "Act. CD3"       = "#fc8d59",
      "Act. CD3+CD28"  = "#d73027"
    )) +
    geom_vline(xintercept = tss_x, linetype = "dashed",
               colour = "black", linewidth = 0.4) +
    annotate("text", x = tss_x, y = Inf, label = "TSS",
             vjust = 1.4, hjust = -0.1, size = 2.2) +
    scale_x_continuous(
      labels = function(x) format(round(x / 1e6, 3), nsmall = 3),
      n.breaks = 4
    ) +
    labs(
      title = sprintf("%s  (FC=%.2f)", gene, fc[gene]),
      x     = sprintf("chr%s (Mb)", as.character(seqnames(prom))),
      y     = "Coverage"
    ) +
    theme_classic(base_size = 9) +
    theme(
      strip.text.y  = element_text(angle = 0, size = 7),
      axis.text.x   = element_text(size = 6, angle = 20, hjust = 1),
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 9)
    )
}

plots <- lapply(sel_idx, plot_promoter)
combined <- Reduce(`+`, plots) + plot_layout(ncol = length(plots))

out_pdf <- "cFOS_ribosomal_promoters.pdf"
pdf(out_pdf, width = 4 * length(plots), height = 6)
print(combined)
dev.off()

message(sprintf("Saved: %s", out_pdf))
