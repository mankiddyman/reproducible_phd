#!/usr/bin/env Rscript
# tpxl6_methyl_tss_check.R -- is x = 0 really the start of each gene in tpxl6_methyl_locus.sh, and how methylated is
# the DNA right at the start?
#
# plot_tpxl6_methyl_locus.R puts x = 0 at the first coding base miniprot found and averages the 2 kb upstream of it.
# That is the promoter only if miniprot found the protein's real start, and even then an active plant promoter is
# unmethylated only in the few hundred bp around the transcription start. For every region (TPXL6 copies and the
# related control gene) this prints and plots:
#   start check  the first query amino acid of the miniprot alignment (1 = the start of the D. regia protein was
#                found) and the BRAKER gene model (results/{sp}/annotation/final/{sp}.final.gff3) that overlaps the
#                miniprot exons most: where its 5' end lies relative to x = 0. A 5' end well upstream of x = 0 means
#                the gene starts earlier and the "2 kb upstream" window lies partly inside it
#   windows      mean CpG methylation, CpGs inside CARP TEs left out, in [-2000,-1000), [-1000,-500), [-500,0) and
#                [0,+500) bp; and the same windows around the BRAKER start when there is one
#   figure       per-CpG methylation from -3 kb to +1.5 kb around x = 0, one row per region; BRAKER 5' end in red
# Uses the BAM choice of plot_tpxl6_methyl_locus.R (methyl_bam_check.tsv: only BAMs whose TEs are methylated).
# FIG_BANNER (environment variable): when set, printed in red at the top right, e.g. to mark a run on test data.
# usage: Rscript tpxl6_methyl_tss_check.R REPO OUT_PREFIX
suppressPackageStartupMessages(library(data.table)); options(width = 220)
a <- commandArgs(TRUE)
if (length(a) < 2) stop("usage: tpxl6_methyl_tss_check.R REPO OUT_PREFIX")
REPO <- a[1]; OUTP <- a[2]
D  <- file.path(REPO, "results/comparative/holocentricity/tpx2/miniprot")
MD <- file.path(REPO, "results/comparative/holocentricity/tpx2/methyl")
Q6 <- "regia_g4881.t1"; Q5 <- "regia_g4900.t1"
X0 <- -3000; X1 <- 1500; SMOOTH <- 100
WIN <- list(`-2000..-1000` = c(-2000, -1000), `-1000..-500` = c(-1000, -500), `-500..0` = c(-500, 0), `0..+500` = c(0, 500))
dir.create(dirname(OUTP), showWarnings = FALSE, recursive = TRUE)
INK <- "#1f1f1e"; INK2 <- "#52514e"; MUTED <- "#8a8984"; RED <- "#b2182b"
GCOL <- c(LTR = "#2a78d6", oTE = "#e87ba4", DNA = "#eda100", rep = "#cfcec9")
alpha <- function(col, a) { m <- col2rgb(col) / 255; rgb(m[1], m[2], m[3], a) }
nice <- function(sp) sub("^Drosera_", "D. ", sp)
hap  <- function(ctg) { h <- regmatches(ctg, regexpr("hap[0-9]+|h[0-9]+tg", ctg)); if (length(h)) sub("^h([0-9]+)tg$", "hap\\1", h) else ctg }
attr1 <- function(x, k) ifelse(grepl(paste0("(^|;)", k, "="), x), sub(paste0(".*(^|;)", k, "=([^;]*).*"), "\\2", x), NA)
gff_read <- function(f, cmd = NULL) {
  if (is.null(cmd)) { if (!file.exists(f) || file.size(f) == 0) return(NULL); x <- readLines(f) } else x <- system(cmd, intern = TRUE)
  x <- x[!startsWith(x, "#") & nzchar(x)]
  if (!length(x)) return(NULL)
  fread(text = x, header = FALSE, sep = "\t", quote = "", col.names = c("contig", "src", "type", "start", "end", "score", "strand", "phase", "attr"))
}
bed_read <- function(f) {
  if (!file.exists(f)) return(NULL)
  b <- fread(cmd = paste("gzip -dc", shQuote(f), "| grep -v '^#'"), header = FALSE, sep = "\t")
  if (!nrow(b)) return(NULL)
  setnames(b, 1:6, c("contig", "begin", "end", "mod", "type", "cov")); b[, pos := begin + 1L][, .(contig, pos, mod, cov)]
}
grp_of <- function(cl) fifelse(startsWith(cl, "Class_I/LTR"), "LTR", fifelse(startsWith(cl, "Class_II"), "DNA",
                      fifelse(startsWith(cl, "Class_I"), "oTE",
                      fifelse(startsWith(cl, "Simple_repeat") | startsWith(cl, "Low_complexity"), NA_character_, "rep"))))
inside <- function(p, s, e) { if (!length(s)) return(rep(FALSE, length(p))); o <- rep(FALSE, length(p)); for (i in seq_along(s)) o <- o | (p >= s[i] & p < e[i]); o }
ovl <- function(a0, a1, b0, b1) sum(pmax(0, pmin(a1, b1) - pmax(a0, b0)))
first_of <- function(x, y) { x <- as.character(x); y <- as.character(y); x[is.na(x)] <- y[is.na(x)]; x }

L <- fread(file.path(D, "tpxl6_miniprot_loci.tsv"))
chk <- if (file.exists(file.path(MD, "methyl_bam_check.tsv"))) fread(file.path(MD, "methyl_bam_check.tsv")) else data.table(species = character(), bam_n = integer(), usable = logical())
sps <- list.dirs(MD, full.names = FALSE, recursive = FALSE); sps <- sps[file.exists(file.path(MD, sps, "regions.tsv"))]
rows <- list(); out <- list()
for (sp in sps) {
  d <- file.path(MD, sp)
  use <- chk[species == sp & usable == TRUE]$bam_n
  bed <- if (length(use) == 2) bed_read(file.path(d, paste0(sp, ".combined.bed.gz"))) else
         if (length(use) == 1) bed_read(file.path(d, sprintf("%s.hap%d.bed.gz", sp, use))) else NULL
  reg <- fread(file.path(d, "regions.tsv"), header = FALSE, col.names = c("set", "contig", "start", "end", "strand"))
  rp <- gff_read(file.path(d, "repeats.gff"))
  if (!is.null(rp)) { rp[, grp := grp_of(first_of(attr1(attr, "classification"), attr1(attr, "Name")))]; rp <- rp[!is.na(grp)] }
  g <- gff_read(file.path(D, paste0(sp, ".gff")))
  bf <- file.path(REPO, "results", sp, "annotation/final", paste0(sp, ".final.gff3"))
  for (i in seq_len(nrow(reg))) {
    r <- reg[i]; flip <- r$strand == "-"; q <- if (r$set == "TPXL6") Q6 else Q5
    tox <- function(p) if (flip) r$end - p else p - r$start
    m <- g[type == "mRNA" & contig == r$contig & start == r$start & end == r$end & startsWith(attr1(attr, "Target"), paste0(q, " "))][1]
    cds <- g[type == "CDS" & attr1(attr, "Parent") == attr1(m$attr, "ID")]
    ex <- data.table(x0 = pmin(tox(cds$start), tox(cds$end)), x1 = pmax(tox(cds$start), tox(cds$end)) + 1,
                     qs = as.integer(tstrsplit(attr1(cds$attr, "Target"), " ")[[2]]))[order(x0)]
    first_aa <- ex$qs[1]
    lrow <- L[species == sp & contig == r$contig & start == r$start & end == r$end & query == q][1]
    # BRAKER model overlapping the miniprot exons most, and where its 5' end lies
    bk <- if (file.exists(bf)) gff_read(cmd = sprintf("awk -F'\\t' '$1 == \"%s\" && $4 <= %d && $5 >= %d' %s", r$contig, r$end + 5000, r$start - 5000, shQuote(bf))) else NULL
    b_id <- NA_character_; b_x5 <- NA_real_; b_str <- NA_character_; b_n <- 0L; b_ncds <- NA_integer_
    if (!is.null(bk)) {
      mr <- bk[type %in% c("mRNA", "transcript")]
      if (nrow(mr)) {
        mr[, id := first_of(attr1(attr, "ID"), attr1(attr, "transcript_id"))]
        bc <- bk[type == "CDS"][, par := attr1(attr, "Parent")]
        mr[, ov := sapply(id, function(z) { cc <- bc[par == z]; if (!nrow(cc)) return(0)
              x0 <- pmin(tox(cc$start), tox(cc$end)); x1 <- pmax(tox(cc$start), tox(cc$end)) + 1
              sum(sapply(seq_along(x0), function(j) ovl(x0[j], x1[j], ex$x0, ex$x1))) })]
        b_n <- nrow(mr)
        best <- mr[which.max(ov)]
        if (best$ov > 0) { b_id <- best$id; b_str <- best$strand; b_x5 <- if (best$strand == "+") tox(best$start) else tox(best$end)
                           b_ncds <- nrow(bc[par == best$id]) }
      }
    }
    te <- if (is.null(rp)) data.table(x0 = numeric(), x1 = numeric()) else rp[contig == r$contig & grp != "rep", .(x0 = pmin(tox(start), tox(end)), x1 = pmax(tox(start), tox(end)) + 1)]
    rr <- if (is.null(rp)) data.table() else rp[contig == r$contig, .(grp, x0 = pmin(tox(start), tox(end)), x1 = pmax(tox(start), tox(end)) + 1)][x1 > X0 & x0 < X1]
    cg <- if (is.null(bed)) data.table(x = numeric(), mod = numeric()) else bed[contig == r$contig, .(x = tox(pos), mod)][x >= -6000 & x < 6000][order(x)]
    cg[, in_te := inside(x, te$x0, te$x1)]
    wmean <- function(w, at = 0) { v <- cg[!in_te & x >= at + w[1] & x < at + w[2]]$mod; c(if (length(v)) round(mean(v), 1) else NA, length(v)) }
    lab <- if (r$set == "TPXL6") paste("TPXL6", hap(r$contig)) else paste("related gene,", hap(r$contig))
    res <- data.table(species = sp, region = lab, contig = r$contig, strand = r$strand,
                      first_query_aa = first_aa, qcov = lrow$qcov, braker_models_near = b_n, braker_best = b_id,
                      braker_cds_exons = b_ncds, braker_5prime_vs_x0 = b_x5)
    for (w in names(WIN)) { v <- wmean(WIN[[w]]); res[[paste0("meth_", w)]] <- v[1]; res[[paste0("n_", w)]] <- v[2] }
    if (!is.na(b_x5) && b_x5 < -100) for (w in names(WIN)) res[[paste0("braker_meth_", w)]] <- wmean(WIN[[w]], b_x5)[1]
    out[[length(out) + 1]] <- res
    if (nrow(cg)) cg[, smooth := sapply(x, function(z) mean(mod[abs(x - z) <= SMOOTH]))]
    rows[[length(rows) + 1]] <- list(sp = sp, lab = lab, ex = ex, rr = rr, cg = cg[x >= X0 & x <= X1], b_x5 = b_x5, has = !is.null(bed))
  }
}
tab <- rbindlist(out, fill = TRUE)
fwrite(tab, paste0(OUTP, ".tsv"), sep = "\t")

n <- length(rows)
draw <- function() {
  par(mai = c(0.7, 2.0, 1.0, 0.3), xpd = NA, family = "sans")
  plot(NA, xlim = c(X0, X1), ylim = c(n + 0.6, 0.45), xaxs = "i", axes = FALSE, xlab = "", ylab = "")
  axis(1, at = seq(X0, X1, 500), labels = seq(X0, X1, 500) / 1000, cex.axis = 0.72, col = "grey55", col.axis = INK2)
  mtext("kb from the first coding base miniprot found (x = 0)", side = 1, line = 2.2, cex = 0.7, col = INK2)
  for (k in seq_len(n)) {
    r <- rows[[k]]; yt <- k - 0.40; yb <- k + 0.12; yg <- k + 0.30; ytr <- function(m) yb - m / 100 * (yb - yt)
    text(X0 - 0.03 * (X1 - X0), k - 0.2, paste(nice(r$sp), r$lab, sep = "\n"), adj = 1, cex = 0.72, col = INK)
    if (nrow(r$rr)) rect(pmax(r$rr$x0, X0), yt, pmin(r$rr$x1, X1), yb, col = sapply(GCOL[r$rr$grp], alpha, 0.28), border = NA)
    segments(X0, c(yt, ytr(50), yb), X1, c(yt, ytr(50), yb), col = c("grey88", "grey80", "grey70"), lty = c(1, 3, 1), lwd = 0.7)
    text(X0 - 0.005 * (X1 - X0), c(yt, yb), c("100%", "0"), adj = c(1, 0.5), cex = 0.5, col = MUTED)
    if (r$has && nrow(r$cg)) { points(r$cg$x, ytr(r$cg$mod), pch = 16, cex = 0.35, col = alpha(INK, 0.5)); lines(r$cg$x, ytr(r$cg$smooth), col = INK, lwd = 1.2) }
    segments(X0, yg, X1, yg, col = "#d6d5d0"); ex <- r$ex[x0 < X1]
    rect(pmax(ex$x0, X0), yg - 0.07, pmin(ex$x1, X1), yg + 0.07, col = INK, border = NA)
    segments(0, yt, 0, yg, lty = 3, col = MUTED)
    if (!is.na(r$b_x5) && r$b_x5 >= X0 && r$b_x5 <= X1) { segments(r$b_x5, yt, r$b_x5, yg + 0.1, col = RED, lwd = 1.3)
      text(r$b_x5, yt - 0.03, "BRAKER start", col = RED, cex = 0.55, adj = c(0.5, 0)) }
  }
  mtext("CpG methylation around the start of each gene", side = 3, line = 3.4, adj = 0, cex = 1.0, font = 2, col = INK)
  mtext(sprintf("Dots: CpGs (%% of reads methylated); line: mean within +-%d bp. Black: coding exons from miniprot; dotted line: x = 0; red: 5' end of the overlapping BRAKER model.", SMOOTH),
        side = 3, line = 2.3, adj = 0, cex = 0.66, col = INK2)
  mtext("Shaded: CARP repeats (blue LTR, pink other retrotransposons, yellow DNA transposons, grey satellite/unclassified).", side = 3, line = 1.5, adj = 0, cex = 0.66, col = INK2)
  if (nzchar(Sys.getenv("FIG_BANNER"))) mtext(Sys.getenv("FIG_BANNER"), side = 3, line = 3.4, adj = 1, cex = 0.95, font = 2, col = RED)
}
dev_png <- function(f, w, h) if (capabilities("cairo")) png(f, width = w, height = h, units = "in", res = 220, type = "cairo") else png(f, width = w, height = h, units = "in", res = 220)
H_IN <- 1.9 + 0.95 * n
dev_png(paste0(OUTP, ".png"), 10.5, H_IN); draw(); invisible(dev.off())
pdf(paste0(OUTP, ".pdf"), 10.5, H_IN); draw(); invisible(dev.off())
cat("== is x = 0 the start of the gene?  (first_query_aa 1 = protein start found; braker_5prime_vs_x0 < 0 = the BRAKER model starts that many bp further upstream)\n")
print(tab[, .(species = nice(species), region, strand, first_query_aa, qcov, braker_models_near, braker_best, braker_cds_exons, braker_5prime_vs_x0)], row.names = FALSE)
cat("\n== mean CpG methylation (%) by window around x = 0, CpGs inside TEs left out (n = CpGs)\n")
wc <- c("species", "region", unlist(lapply(names(WIN), function(w) c(paste0("meth_", w), paste0("n_", w)))))
pr <- tab[, ..wc]; pr$species <- nice(pr$species); print(pr, row.names = FALSE)
bc <- grep("^braker_meth_", names(tab), value = TRUE)
if (length(bc)) { cat("\n== the same windows around the BRAKER start, where it lies > 100 bp upstream of x = 0\n")
  pb <- tab[!is.na(get(bc[1])), c("species", "region", "braker_5prime_vs_x0", bc), with = FALSE]; pb$species <- nice(pb$species); print(pb, row.names = FALSE) }
cat(sprintf("wrote %s.png / .pdf / .tsv\n", OUTP))
