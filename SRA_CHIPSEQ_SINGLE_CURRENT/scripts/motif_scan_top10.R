suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(rtracklayer)
})

BASE <- "/Volumes/work/obst/GSE116695_new"

# ── Rebuild top 10 peaks ──────────────────────────────────────────────────────
gtf <- import(file.path(BASE, "Homo_sapiens.GRCh38.115.gtf.gz"),
              format = "gtf", feature.type = "gene",
              colnames = c("gene_name", "gene_biotype"))
ribo_gr <- gtf[grepl("^RP[LS][0-9]", gtf$gene_name) &
               gtf$gene_biotype == "protein_coding" &
               as.character(seqnames(gtf)) %in% as.character(1:22)]
ribo_gr <- ribo_gr[!duplicated(ribo_gr$gene_name)]
exclude  <- c("RPS6KA2","RPS6KA5","RPS6KB1","RPS6KB2","RPS6KC1","RPS19BP1")
ribo_gr  <- ribo_gr[!ribo_gr$gene_name %in% exclude]

tss_gr <- GRanges(
  seqnames = seqnames(ribo_gr),
  ranges   = IRanges(
    start = ifelse(as.character(strand(ribo_gr)) == "+", start(ribo_gr), end(ribo_gr)),
    width = 1),
  strand = strand(ribo_gr), gene = ribo_gr$gene_name)
promoter_gr <- suppressWarnings(promoters(tss_gr, upstream = 2000, downstream = 200))

np_files <- list.files(file.path(BASE, "peak_calls"),
                       pattern = "_peaks\\.narrowPeak$", full.names = TRUE)
peaks_all <- do.call(c, lapply(np_files, function(f) {
  pk <- import(f, format = "narrowPeak")
  pk$comparison <- sub("_peaks\\.narrowPeak$", "", basename(f))
  pk
}))
ol       <- findOverlaps(peaks_all, promoter_gr)
hit      <- peaks_all[queryHits(ol)]
hit$gene <- tss_gr$gene[subjectHits(ol)]
hit      <- hit[!duplicated(paste(hit$comparison, seqnames(hit), start(hit)))]
hit      <- hit[order(-hit$signalValue)]
top10    <- head(hit, 10)
top10$tf <- sub("__.*", "", top10$comparison)

# ── JASPAR PFMs as Biostrings PWMatrix objects ────────────────────────────────
# NFATC2 / NFAT1 — MA0152.1 (consensus fwd strand: TTTTCCA; RC = TGGAAAA)
pfm_nfat1 <- rbind(
  A = c(  3,  1,  1,  0,  1,  0, 18),
  C = c(  1,  2,  1,  0, 25, 26,  3),
  G = c(  2,  2,  0,  0,  0,  0,  1),
  T = c( 20, 21, 24, 26,  0,  0,  4)
)

# JUNB — MA0490.3 (ChIP-seq; AP-1 consensus ATGACTCAT)
pfm_junb <- rbind(
  A = c(56607,   538,   770, 78101,  4755,  1999,  5243, 77675,  3722),
  C = c( 7817,   654,   550,   523, 59784,   289, 73199,   588,  8876),
  G = c(11566,   322, 74757,   550, 10519,   913,   617,   900,  7508),
  T = c( 3866, 78342,  3779,   682,  4798, 76655,   797,   693, 59750)
)

# Convert PFM -> PWM (plain matrix for Biostrings::matchPWM): log-odds with pseudocount
pfm_to_pwm <- function(pfm, pseudo = 0.8) {
  pfm <- pfm + pseudo
  ppm <- sweep(pfm, 2, colSums(pfm), "/")
  log2(ppm / 0.25)   # log-odds vs uniform background; rownames A,C,G,T preserved
}

pwm_nfat1 <- pfm_to_pwm(pfm_nfat1)
pwm_junb  <- pfm_to_pwm(pfm_junb)
pwms <- list(NFAT1 = pwm_nfat1, JUNB = pwm_junb)

genome <- BSgenome.Hsapiens.UCSC.hg38

# ── Helper: scan one peak with its TF's PWM ───────────────────────────────────
scan_peak <- function(gr_peak, tf, min_score = "80%") {
  pwm <- pwms[[tf]]
  if (is.null(pwm)) return(list(found = NA, n = NA, best = NA))

  # add chr prefix for UCSC genome
  chrom <- paste0("chr", as.character(seqnames(gr_peak)))
  if (!chrom %in% seqnames(genome)) return(list(found = FALSE, n = 0, best = NA))

  seq_fwd <- getSeq(genome,
                    GRanges(chrom, ranges(gr_peak), strand = "+"))[[1]]
  seq_rev <- reverseComplement(seq_fwd)

  hits_fwd <- matchPWM(pwm, seq_fwd, min.score = min_score, with.score = TRUE)
  hits_rev <- matchPWM(pwm, seq_rev, min.score = min_score, with.score = TRUE)

  n <- length(hits_fwd) + length(hits_rev)
  if (n == 0) return(list(found = FALSE, n = 0, best = NA))

  all_scores <- c(mcols(hits_fwd)$score, mcols(hits_rev)$score)
  best_score <- max(all_scores)
  max_possible <- sum(apply(pwm, 2, max))
  list(found = TRUE, n = n, best = round(best_score / max_possible * 100, 1))
}

# ── Scan all 10 peaks ─────────────────────────────────────────────────────────
cat(sprintf("Scanning top 10 peaks with JASPAR PWMs (NFAT1: MA0152.1, JUNB: MA0490.3)\n"))
cat(sprintf("Threshold: 80%% of maximum possible PWM score (both strands)\n\n"))

rows <- lapply(seq_along(top10), function(i) {
  tf   <- as.character(top10$tf[i])
  gene <- as.character(top10$gene[i])
  comp <- as.character(top10$comparison[i])
  pw   <- as.integer(end(top10[i]) - start(top10[i]) + 1)
  res  <- scan_peak(top10[i], tf, min_score = "80%")

  data.frame(
    rank        = i,
    gene        = gene,
    tf          = tf,
    peak_width  = pw,
    motif_found = res$found,
    n_hits      = res$n,
    best_pct_max = if (is.na(res$best)) NA_character_ else paste0(res$best, "%"),
    stringsAsFactors = FALSE
  )
})

res_df <- do.call(rbind, rows)
print(res_df, row.names = FALSE)

# ── Build per-hit table with genomic coordinates ─────────────────────────────
hit_rows <- list()

for (i in seq_along(top10)) {
  tf    <- as.character(top10$tf[i])
  gene  <- as.character(top10$gene[i])
  comp  <- as.character(top10$comparison[i])
  chrom_ens <- as.character(seqnames(top10[i]))
  chrom_ucsc <- paste0("chr", chrom_ens)
  pk_start  <- start(top10[i])   # 1-based genomic
  pk_end    <- end(top10[i])
  L         <- pk_end - pk_start + 1
  sig       <- round(top10$signalValue[i], 3)
  pwm       <- pwms[[tf]]

  if (is.null(pwm) || !chrom_ucsc %in% seqnames(genome)) next

  seq_fwd <- getSeq(genome, GRanges(chrom_ucsc, IRanges(pk_start, pk_end), strand = "+"))[[1]]
  seq_rev <- reverseComplement(seq_fwd)
  m_len   <- ncol(pwm)

  hits_f <- matchPWM(pwm, seq_fwd, min.score = "80%", with.score = TRUE)
  hits_r <- matchPWM(pwm, seq_rev, min.score = "80%", with.score = TRUE)

  # Forward hits: genomic_start = pk_start + (within_peak_pos - 1)
  if (length(hits_f) > 0) {
    for (j in seq_along(hits_f)) {
      p <- start(hits_f)[j]
      hit_rows[[length(hit_rows) + 1]] <- data.frame(
        peak_rank        = i,
        gene             = gene,
        tf               = tf,
        comparison       = comp,
        peak_signalValue = sig,
        chrom            = chrom_ens,
        motif_start      = pk_start + p - 1,
        motif_end        = pk_start + p - 1 + m_len - 1,
        strand           = "+",
        sequence         = as.character(hits_f)[j],
        score            = round(mcols(hits_f)$score[j], 3),
        pct_max_score    = round(mcols(hits_f)$score[j] / sum(apply(pwm, 2, max)) * 100, 1),
        stringsAsFactors = FALSE
      )
    }
  }

  # Reverse hits: hit at pos p in RC → genomic_start = pk_start + (L - p - m_len + 1) - 1 + 1
  #   = pk_start + L - p - m_len + 1
  if (length(hits_r) > 0) {
    for (j in seq_along(hits_r)) {
      p <- start(hits_r)[j]
      g_start <- pk_start + L - p - m_len + 1
      g_end   <- pk_start + L - p
      hit_rows[[length(hit_rows) + 1]] <- data.frame(
        peak_rank        = i,
        gene             = gene,
        tf               = tf,
        comparison       = comp,
        peak_signalValue = sig,
        chrom            = chrom_ens,
        motif_start      = g_start,
        motif_end        = g_end,
        strand           = "-",
        sequence         = as.character(hits_r)[j],
        score            = round(mcols(hits_r)$score[j], 3),
        pct_max_score    = round(mcols(hits_r)$score[j] / sum(apply(pwm, 2, max)) * 100, 1),
        stringsAsFactors = FALSE
      )
    }
  }
}

hits_df <- do.call(rbind, hit_rows)
hits_df <- hits_df[order(hits_df$peak_rank, hits_df$motif_start), ]

cat("\nMotif hits with genomic coordinates:\n\n")
print(hits_df, row.names = FALSE)

# Save TSV
out_file <- file.path(BASE, "top_motif_hits.tsv")
write.table(hits_df, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nSaved %d rows to %s\n", nrow(hits_df), out_file))
