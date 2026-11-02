#!/usr/bin/env Rscript
# plot_tpxl6_methyl_locus.R -- CpG methylation over every TPXL6 copy and an expressed relative (tpxl6_methyl_locus.sh)
#
# Question: is D. binata TPXL6 silenced by DNA methylation (methylated promoter, like a TE), or simply not switched on
# in leaves (unmethylated promoter)? Each species' expressed TPXL1/5-group gene is the "unmethylated promoter"
# reference, and the TEs inside the loci are the "methylated" reference.
# Per BAM check first: TEs are methylated in any native plant DNA, so a BAM whose reads show TEs below MIN_TE % is
# reported as not usable (typically an amplified library: PCR copies carry no methylation) and left out of the plot.
# One row per region, 5'->3' along the protein, x = 0 at the first coding base miniprot found (exon 2 for a copy
# whose first exon was not found), UP bp upstream to DOWN bp past the last coding base:
#   dots   % of reads methylated at each CpG (pb-CpG-tools count mode, >= 3 reads); line: mean over +-200 bp
#   shade  CARP repeats (colours as in plot_tpxl6_te_loci.R); grey band = the 2 kb upstream window
#   below  coding exons (black) on the locus line
# Columns: mean CpG methylation in the 2 kb upstream, in coding exons, in intron sequence outside TEs, in TEs; median
# reads per CpG; the species' leaf-RNA spliced reads (tpxl6_rna_kmer.sh). Also writes methyl_summary.tsv.
# usage: Rscript plot_tpxl6_methyl_locus.R MINIPROT_DIR METHYL_DIR SPECIES_CSV OUT_PREFIX [RNA_KMER_SUMMARY_TSV]
suppressPackageStartupMessages(library(data.table)); options(width = 220)
a <- commandArgs(TRUE)
if (length(a) < 4) stop("usage: plot_tpxl6_methyl_locus.R MINIPROT_DIR METHYL_DIR SPECIES_CSV OUT_PREFIX [RNA_TSV]")
D <- a[1]; MD <- a[2]; SPCSV <- a[3]; OUTP <- a[4]; RNA <- if (length(a) >= 5) a[5] else ""
Q6 <- "regia_g4881.t1"; Q5 <- "regia_g4900.t1"
UP <- 5000; DOWN <- 2000; UPWIN <- 2000; SMOOTH <- 200; MINEX <- 250; MIN_TE <- 40
NAME_IN <- 2.0; TRACK_IN <- 6.2; COLS_IN <- c(0.8, 0.7, 0.85, 0.65, 0.75, 1.25)
dir.create(dirname(OUTP), showWarnings = FALSE, recursive = TRUE)

ORDER <- c("Drosera_regia", "Drosera_paradoxa", "Drosera_roseana", "Drosera_scorpioides", "Drosera_binata")
cen <- c(Drosera_regia = "holocentric", Drosera_capensis = "monocentric")
if (file.exists(SPCSV)) { s <- fread(SPCSV, select = c("species_id", "centromere")); cen[s$species_id] <- s$centromere }
CENCOL <- c(holocentric = "#1b9e77", monocentric = "#d95f02")
GRP  <- c(LTR = "LTR retrotransposons", oTE = "LINE, SINE, other TEs", DNA = "DNA transposons", rep = "satellite / unclassified repeat")
GCOL <- c(LTR = "#2a78d6", oTE = "#e87ba4", DNA = "#eda100", rep = "#cfcec9")
INK <- "#1f1f1e"; INK2 <- "#52514e"; MUTED <- "#8a8984"; FLANK <- "#d6d5d0"; UPBAND <- "#ecebe6"
alpha <- function(col, a) { m <- col2rgb(col) / 255; rgb(m[1], m[2], m[3], a) }
nice <- function(sp) sub("^Drosera_", "D. ", sp)
hap  <- function(ctg) { h <- regmatches(ctg, regexpr("hap[0-9]+|h[0-9]+tg", ctg)); if (length(h)) sub("^h([0-9]+)tg$", "hap\\1", h) else ctg }
attr1 <- function(x, k) ifelse(grepl(paste0("(^|;)", k, "="), x), sub(paste0(".*(^|;)", k, "=([^;]*).*"), "\\2", x), NA)
gff_read <- function(f) {
  if (!file.exists(f) || file.size(f) == 0) return(NULL)
  x <- readLines(f); x <- x[!startsWith(x, "#") & nzchar(x)]
  if (!length(x)) return(NULL)
  fread(text = x, header = FALSE, sep = "\t", quote = "",
        col.names = c("contig", "src", "type", "start", "end", "score", "strand", "phase", "attr"))
}
bed_read <- function(f) {                       # pb-CpG-tools bed.gz: first 6 columns are shared by both pileup modes
  if (!file.exists(f)) return(NULL)
  b <- fread(cmd = paste("gzip -dc", shQuote(f), "| grep -v '^#'"), header = FALSE, sep = "\t")
  if (!nrow(b)) return(NULL)
  setnames(b, 1:6, c("contig", "begin", "end", "mod", "type", "cov"))
  b[, pos := begin + 1L][, .(contig, pos, mod, cov)]
}
grp_of <- function(cl) fifelse(startsWith(cl, "Class_I/LTR"), "LTR", fifelse(startsWith(cl, "Class_II"), "DNA",
                      fifelse(startsWith(cl, "Class_I"), "oTE",
                      fifelse(startsWith(cl, "Simple_repeat") | startsWith(cl, "Low_complexity"), NA_character_, "rep"))))
inside <- function(p, s, e) { if (!length(s)) return(rep(FALSE, length(p))); o <- rep(FALSE, length(p)); for (i in seq_along(s)) o <- o | (p >= s[i] & p < e[i]); o }
mn <- function(v) if (length(v)) mean(v) else NA_real_

L <- fread(file.path(D, "tpxl6_miniprot_loci.tsv"))
rna <- if (nzchar(RNA) && file.exists(RNA)) fread(RNA) else NULL
sps <- list.dirs(MD, full.names = FALSE, recursive = FALSE)
sps <- sps[file.exists(file.path(MD, sps, "regions.tsv"))]
sps <- c(intersect(ORDER, sps), setdiff(sps, ORDER))
if (!length(sps)) stop("no species folders with regions.tsv in ", MD)

rows <- list(); checks <- list(); notes <- character()
for (sp in sps) {
  d <- file.path(MD, sp)
  reg <- fread(file.path(d, "regions.tsv"), header = FALSE, col.names = c("set", "contig", "start", "end", "strand"))
  rp <- gff_read(file.path(d, "repeats.gff"))
  if (!is.null(rp)) { rp[, grp := grp_of(fifelse(is.na(attr1(attr, "classification")), attr1(attr, "Name"), attr1(attr, "classification")))]; rp <- rp[!is.na(grp)] }
  bams <- if (file.exists(file.path(d, "bams.tsv"))) fread(file.path(d, "bams.tsv")) else data.table(n = integer(), bam = character())
  # per-BAM check: methylation of the TEs in these regions
  te_iv <- if (is.null(rp)) data.table() else rp[grp != "rep", .(contig, s = start, e = end + 1)]
  te_meth <- function(b) { if (is.null(b) || !nrow(te_iv)) return(NA_real_)
    v <- unlist(lapply(split(te_iv, te_iv$contig), function(t) { x <- b[contig == t$contig[1]]; x$mod[inside(x$pos, t$s, t$e)] }))
    if (length(v) < 20) NA_real_ else mean(v) }
  use <- c()
  for (h in 1:2) {
    b <- bed_read(file.path(d, sprintf("%s.hap%d.bed.gz", sp, h)))
    if (is.null(b)) next
    tm <- te_meth(b); ok <- !is.na(tm) && tm >= MIN_TE
    nm <- if (nrow(bams[n == h])) basename(bams[n == h]$bam) else sprintf("BAM %d", h)
    checks[[length(checks) + 1]] <- data.table(species = sp, bam_n = h, bam = nm, CpGs = nrow(b),
                                               TE_CpG_methylation = round(tm, 1), usable = ok)
    if (ok) use <- c(use, h)
  }
  bed <- if (length(use) == 2) bed_read(file.path(d, paste0(sp, ".combined.bed.gz"))) else
         if (length(use) == 1) bed_read(file.path(d, sprintf("%s.hap%d.bed.gz", sp, use))) else NULL
  if (is.null(bed)) notes <- c(notes, sprintf("%s: no BAM with methylated TEs, so no usable methylation", nice(sp)))
  g <- gff_read(file.path(D, paste0(sp, ".gff")))
  for (i in seq_len(nrow(reg))) {
    r <- reg[i]; flip <- r$strand == "-"; q <- if (r$set == "TPXL6") Q6 else Q5
    tox <- function(p) if (flip) r$end - p else p - r$start                      # 1-based position -> x
    m <- g[type == "mRNA" & contig == r$contig & start == r$start & end == r$end & startsWith(attr1(attr, "Target"), paste0(q, " "))][1]
    cds <- g[type == "CDS" & attr1(attr, "Parent") == attr1(m$attr, "ID")]
    ex <- data.table(x0 = pmin(tox(cds$start), tox(cds$end)), x1 = pmax(tox(cds$start), tox(cds$end)) + 1)[order(x0)]
    span <- max(ex$x1)
    intr <- if (nrow(ex) > 1) data.table(a = ex$x1[-nrow(ex)], b = ex$x0[-1]) else data.table(a = numeric(), b = numeric())
    rr <- if (is.null(rp)) data.table() else rp[contig == r$contig]
    if (nrow(rr)) { rr[, `:=`(x0 = pmin(tox(start), tox(end)), x1 = pmax(tox(start), tox(end)) + 1)]; rr <- rr[x1 > -UP & x0 < span + DOWN] }
    cg <- if (is.null(bed)) data.table() else bed[contig == r$contig, .(x = tox(pos), mod, cov)][x >= -UP & x < span + DOWN][order(x)]
    te <- if (nrow(rr)) rr[grp != "rep"] else data.table(x0 = numeric(), x1 = numeric())
    if (nrow(cg)) {
      cg[, `:=`(in_te = inside(x, te$x0, te$x1), in_ex = inside(x, ex$x0, ex$x1), in_intr = inside(x, intr$a, intr$b))]
      cg[, smooth := sapply(x, function(z) mean(mod[abs(x - z) <= SMOOTH]))]
    }
    w <- function(f) if (!nrow(cg)) NA_real_ else mn(cg[f(cg)]$mod)
    lab <- if (r$set == "TPXL6") hap(r$contig) else paste0("related gene, ", hap(r$contig))
    rows[[length(rows) + 1]] <- list(sp = sp, set = r$set, lab = lab, contig = r$contig, ex = ex, span = span, rr = rr, cg = cg,
      has = !is.null(bed),
      up = w(function(t) t$x >= -UPWIN & t$x < 0 & !t$in_te), exon = w(function(t) t$in_ex),
      intron = w(function(t) t$in_intr & !t$in_te), tes = w(function(t) t$in_te),
      n_up = if (nrow(cg)) sum(cg$x >= -UPWIN & cg$x < 0 & !cg$in_te) else 0L,
      cov = if (nrow(cg)) median(cg$cov) else NA_real_)
  }
}
chk <- rbindlist(checks)
summ <- rbindlist(lapply(rows, function(r) data.table(species = r$sp, region = r$set, copy = r$lab, contig = r$contig,
  upstream2kb = round(r$up, 1), upstream_CpGs = r$n_up, exons = round(r$exon, 1), introns_nonTE = round(r$intron, 1),
  TEs = round(r$tes, 1), median_reads_per_CpG = r$cov)))
fwrite(summ, file.path(MD, "methyl_summary.tsv"), sep = "\t")
if (nrow(chk)) fwrite(chk, file.path(MD, "methyl_bam_check.tsv"), sep = "\t")

rna_of <- function(sp, which) {
  if (is.null(rna)) return("")
  w <- if (which == "TPXL6") "TPXL6" else "control"
  x <- rna[species == sp & set == w]
  if (!nrow(x)) "no RNA" else format(x$junction_reads[1], big.mark = ",")
}
pct <- function(v) if (is.na(v)) "-" else sprintf("%.0f%%", v)

n <- length(rows); maxspan <- max(sapply(rows, `[[`, "span"))
X0 <- -UP; X1 <- ceiling((maxspan + DOWN) / 5000) * 5000
u <- (X1 - X0) / TRACK_IN
xc <- X1 + (0.35 + cumsum(COLS_IN)) * u
W_IN <- NAME_IN + TRACK_IN + 0.35 + sum(COLS_IN) + 0.2
NL <- length(notes); MAI_TOP <- (6.0 + NL) * 0.2
H_IN <- MAI_TOP + 0.95 + 0.95 * n
draw <- function() {
  par(mai = c(0.75, 0.1, MAI_TOP, 0.1), xpd = NA, family = "sans")
  plot(NA, xlim = c(X0 - NAME_IN * u, X1 + (0.35 + sum(COLS_IN) + 0.1) * u), ylim = c(n + 0.6, 0.45), xaxs = "i",
       axes = FALSE, xlab = "", ylab = "")
  at <- c(X0 / 1000, pretty(c(0, X1 / 1000))); at <- unique(at[at * 1000 >= X0 & at * 1000 <= X1])
  axis(1, at = at * 1000, labels = at, cex.axis = 0.72, col = "grey55", col.axis = INK2, line = 0.2)
  mtext("kb from the first coding base (negative = upstream)", side = 1, line = 2.3, at = (X0 + X1) / 2, cex = 0.7, col = INK2)
  yl <- function(L) { u4 <- par("usr"); u4[4] - L * par("csi") * (u4[3] - u4[4]) / par("pin")[2] }   # L lines above the plot
  hdr <- c("2 kb\nupstream", "coding\nexons", "introns,\nnot TE", "TEs", "reads\nper CpG", "leaf RNA, spliced\nreads (all copies)")
  text(xc, yl(1.0), hdr, adj = c(1, 0.5), cex = 0.66, font = 2, col = INK)
  last <- ""
  for (k in seq_len(n)) {
    r <- rows[[k]]; yt <- k - 0.40; yb <- k + 0.12; yg <- k + 0.30
    ytr <- function(m) yb - m / 100 * (yb - yt)
    if (r$sp != last) {
      if (k > 1) segments(X0 - (NAME_IN - 0.05) * u, k - 0.5, xc[6], k - 0.5, col = "grey85")
      cc <- if (is.na(cen[r$sp])) NA else CENCOL[cen[r$sp]]
      if (!is.na(cc)) points(X0 - (NAME_IN - 0.08) * u, k - 0.32, pch = 15, cex = 0.95, col = cc)
      text(X0 - (NAME_IN - 0.18) * u, k - 0.32, nice(r$sp), adj = 0, font = 3, cex = 0.82, col = INK)
      if (!is.na(cen[r$sp])) text(X0 - (NAME_IN - 0.18) * u, k - 0.08, cen[r$sp], adj = 0, cex = 0.66, col = INK2)
      last <- r$sp
    }
    text(X0 - (NAME_IN - 0.18) * u, k + 0.18, r$lab, adj = 0, cex = 0.7, font = 2, col = INK)
    rect(-UPWIN, yt, 0, yb, col = UPBAND, border = NA)                                  # 2 kb upstream window
    if (nrow(r$rr)) { rr <- r$rr[order(grp != "rep")]
      rect(pmax(rr$x0, -UP), yt, pmin(rr$x1, r$span + DOWN), yb, col = sapply(GCOL[rr$grp], alpha, 0.28), border = NA) }
    segments(-UP, c(yt, ytr(50), yb), r$span + DOWN, c(yt, ytr(50), yb), col = c("grey88", "grey80", "grey70"),
             lty = c(1, 3, 1), lwd = 0.7)
    text(-UP - 0.05 * u, c(yt, yb), c("100%", "0"), adj = c(1, 0.5), cex = 0.5, col = MUTED)
    if (r$has && nrow(r$cg)) {
      points(r$cg$x, ytr(r$cg$mod), pch = 16, cex = 0.28, col = alpha(INK, 0.45))
      lines(r$cg$x, ytr(r$cg$smooth), col = INK, lwd = 1.2)
    } else text((X0 + X1) / 2, (yt + yb) / 2, "no usable methylation calls", cex = 0.7, font = 3, col = MUTED)
    segments(-UP, yg, r$span + DOWN, yg, col = FLANK, lwd = 1); segments(0, yg, r$span, yg, col = MUTED, lwd = 1.1)
    rect(r$ex$x0, yg - 0.07, pmax(r$ex$x1, r$ex$x0 + MINEX), yg + 0.07, col = INK, border = NA)
    vals <- c(pct(r$up), pct(r$exon), pct(r$intron), pct(r$tes), if (is.na(r$cov)) "-" else sprintf("%.0f", r$cov))
    text(xc[1:5], k - 0.14, vals, adj = c(1, 0.5), cex = c(0.78, 0.72, 0.72, 0.72, 0.72), font = c(2, 1, 1, 1, 1), col = INK)
    ks <- which(sapply(rows, function(z) z$sp == r$sp && z$set == r$set))         # RNA: one count per species and gene
    if (k == min(ks)) {
      yc <- mean(range(ks)) - 0.14
      if (length(ks) > 1) { xb <- xc[6] - 0.62 * u; yb <- c(min(ks) - 0.45, max(ks) + 0.17)
        segments(xb, yb[1], xb, yb[2], col = "grey65", lwd = 0.9); segments(xb, yb, xb - 0.06 * u, yb, col = "grey65", lwd = 0.9) }
      text(xc[6], yc, rna_of(r$sp, r$set), adj = c(1, 0.5), cex = 0.72, col = INK)
    }
  }
  xl <- X0 - (NAME_IN - 0.18) * u
  mtext("CpG methylation over TPXL6 (PacBio HiFi 5mC)", side = 3, line = 4.75 + NL, at = xl, adj = 0, cex = 1.05, font = 2, col = INK)
  mtext(sprintf("Each dot is one CpG (%% of reads methylated, >= 3 reads); line = mean within +-%d bp. Shaded: CARP repeats; grey band: the 2 kb upstream window. Black: coding exons.", SMOOTH),
        side = 3, line = 3.65 + NL, at = xl, adj = 0, cex = 0.68, col = INK2)
  mtext("x = 0 is the first coding base found (exon 2 for a copy whose first exon was not found). Upstream % excludes CpGs inside TEs; BAMs whose TEs are unmethylated are left out.",
        side = 3, line = 2.8 + NL, at = xl, adj = 0, cex = 0.68, col = INK2)
  for (j in seq_along(notes)) mtext(notes[j], side = 3, line = 2.8 + NL - j, at = xl, adj = 0, cex = 0.68, col = "#b2182b")
  leg <- function(y, labs, cols) { x <- X0
    for (i in seq_along(labs)) { points(x, y, pch = 15, cex = 1.5, col = cols[i]); text(x + 0.11 * u, y, labs[i], adj = c(0, 0.5), cex = 0.68, col = INK)
      x <- x + 0.11 * u + strwidth(labs[i], cex = 0.68) + 0.32 * u } }
  leg(yl(1.85), GRP[c("LTR", "oTE", "DNA")], sapply(GCOL[c("LTR", "oTE", "DNA")], alpha, 0.5))
  leg(yl(0.95), c(GRP[["rep"]], "2 kb upstream window", "coding exon"), c(alpha(GCOL[["rep"]], 0.5), UPBAND, INK))
}
dev_png <- function(f, w, h) if (capabilities("cairo")) png(f, width = w, height = h, units = "in", res = 220, type = "cairo") else png(f, width = w, height = h, units = "in", res = 220)
dev_png(paste0(OUTP, ".png"), W_IN, H_IN); draw(); invisible(dev.off())
pdf(paste0(OUTP, ".pdf"), W_IN, H_IN); draw(); invisible(dev.off())
cat("== per-BAM check: CpG methylation of the TEs in these regions (native plant DNA: high; amplified library: near 0)\n")
if (nrow(chk)) print(chk, row.names = FALSE) else cat("   no per-BAM pileups found\n")
cat("\n== mean CpG methylation (%) per region\n"); print(summ, row.names = FALSE)
if (length(notes)) cat("\n!!", paste(notes, collapse = "\n!! "), "\n")
cat(sprintf("wrote %s.png / .pdf (%.1f x %.1f in), %s and methyl_bam_check.tsv\n", OUTP, W_IN, H_IN, file.path(MD, "methyl_summary.tsv")))
