#!/usr/bin/env Rscript
# plot_tpxl6_te_loci.R -- transposable elements in and around every TPXL6 copy (CARP annotation)
#
# Tests whether the TE-filled introns of D. binata TPXL6 could explain why it is silent in leaves: do the expressed
# holocentric copies carry TEs in their introns too, and does binata alone have a TE just upstream of the gene?
# One row per TPXL6 copy (assigned TPXL6 in tpxl6_miniprot_loci.tsv), 5'->3' along the protein, x = 0 at the
# first coding base miniprot found, from UP bp upstream to DOWN bp past the last coding base. Repeats are the
# CARP Level-1 features kept by tpxl6_te_loci.sh (one classification per base), grouped:
#   LTR retrotransposons    Class_I/LTR/...
#   LINE, SINE, other TEs   every other Class_I path (LINE, SINE, DIRS, Penelope, ...) and Class_I/Class_II
#   DNA transposons         Class_II/... (TIR, MITE, Helitron)
#   grey                    Satellite, rDNA, Unknown: drawn, but not counted as TE
#   Simple_repeat and Low_complexity are left out.
# Columns: summed intron length; % of intron bp in TEs; % of the UPWIN bp upstream of the first coding base in TEs;
# the species' leaf-RNA reads spanning a TPXL6 splice junction (all copies together, from tpxl6_rna_kmer.sh).
# Also writes te_loci_summary.tsv (one line per copy) into TE_DIR.
# usage: Rscript plot_tpxl6_te_loci.R MINIPROT_DIR TE_DIR SPECIES_CSV OUT_PREFIX [RNA_KMER_SUMMARY_TSV]
suppressPackageStartupMessages(library(data.table)); options(width = 220)
a <- commandArgs(TRUE)
if (length(a) < 4) stop("usage: plot_tpxl6_te_loci.R MINIPROT_DIR TE_DIR SPECIES_CSV OUT_PREFIX [RNA_TSV]")
D <- a[1]; TED <- a[2]; SPCSV <- a[3]; OUTP <- a[4]; RNA <- if (length(a) >= 5) a[5] else ""
Q6 <- "regia_g4881.t1"
UP <- 5000; DOWN <- 2000; UPWIN <- 2000; MINEX <- 250    # bp
NAME_IN <- 1.9; TRACK_IN <- 6.2; COLS_IN <- c(0.85, 0.95, 1.05, 1.25)   # inches: names, track, 4 columns
dir.create(dirname(OUTP), showWarnings = FALSE, recursive = TRUE)

ORDER <- c("Nepenthes_gracilis", "Dionaea_muscipula", "Drosera_regia", "Drosera_paradoxa", "Drosera_roseana",
           "Drosera_scorpioides", "Drosera_binata")
NICE <- c(Nepenthes_gracilis = "Nepenthes gracilis", Dionaea_muscipula = "Dionaea muscipula")
cen <- c(Drosera_regia = "holocentric", Drosera_capensis = "monocentric")
if (file.exists(SPCSV)) { s <- fread(SPCSV, select = c("species_id", "centromere")); cen[s$species_id] <- s$centromere }
CENCOL <- c(holocentric = "#1b9e77", monocentric = "#d95f02")          # as in the gene-structure figure
# TE groups: three hues that stay apart for every pair, with colour-vision deficiency too (dataviz validator,
# all pairs: CVD dE 13.0, normal 19.6), and none of them the centromere colours; grey = repeat, not a TE
GRP  <- c(LTR = "LTR retrotransposons", oTE = "LINE, SINE, other TEs", DNA = "DNA transposons",
          rep = "satellite / unclassified repeat")
GCOL <- c(LTR = "#2a78d6", oTE = "#e87ba4", DNA = "#eda100", rep = "#cfcec9")
INK <- "#1f1f1e"; INK2 <- "#52514e"; MUTED <- "#8a8984"; LINE <- "#8a8984"; FLANK <- "#d6d5d0"

nice <- function(sp) if (sp %in% names(NICE)) NICE[[sp]] else sub("^Drosera_", "D. ", sp)
hap  <- function(ctg) { h <- regmatches(ctg, regexpr("hap[0-9]+|h[0-9]+tg", ctg)); if (length(h)) sub("^h([0-9]+)tg$", "hap\\1", h) else "" }
attr1 <- function(x, k) ifelse(grepl(paste0("(^|;)", k, "="), x), sub(paste0(".*(^|;)", k, "=([^;]*).*"), "\\2", x), NA)
merge_iv <- function(s, e) {                    # union of half-open intervals [s, e)
  if (!length(s)) return(list(s = numeric(), e = numeric()))
  o <- order(s); s <- s[o]; e <- e[o]; rs <- s[1]; re <- e[1]
  for (i in seq_along(s)[-1]) {
    if (s[i] <= re[length(re)]) re[length(re)] <- max(re[length(re)], e[i]) else { rs <- c(rs, s[i]); re <- c(re, e[i]) }
  }
  list(s = rs, e = re)
}
covered <- function(u, a, b) if (!length(u$s)) 0 else sum(pmax(0, pmin(u$e, b) - pmax(u$s, a)))
gff_read <- function(f) {
  x <- readLines(f); x <- x[!startsWith(x, "#") & nzchar(x)]
  if (!length(x)) return(NULL)
  fread(text = x, header = FALSE, sep = "\t", quote = "",
        col.names = c("contig", "src", "type", "start", "end", "score", "strand", "phase", "attr"))
}

L  <- fread(file.path(D, "tpxl6_miniprot_loci.tsv"))
cp_all <- L[query == Q6 & assigned == "TPXL6"]
st <- if (file.exists(file.path(TED, "status.tsv"))) fread(file.path(TED, "status.tsv"), sep = "\t") else
        data.table(species = character(), carp = character(), note = character())
rna <- if (nzchar(RNA) && file.exists(RNA)) fread(RNA) else NULL
sp_all <- c(intersect(ORDER, unique(cp_all$species)), setdiff(unique(cp_all$species), ORDER))

te_of <- function(sp) {                         # CARP features of a species, grouped; NULL = no usable CARP run
  if (!nrow(st[species == sp & carp == "yes"])) return(NULL)
  f <- file.path(TED, paste0(sp, ".te.gff"))
  g <- if (file.exists(f) && file.size(f) > 0) gff_read(f) else NULL
  if (is.null(g)) return(data.table(contig = character(), start = integer(), end = integer(), grp = character()))
  cl <- attr1(g$attr, "classification"); nm <- attr1(g$attr, "Name"); cl[is.na(cl)] <- nm[is.na(cl)]
  g[, grp := fifelse(startsWith(cl, "Class_I/LTR"), "LTR",
             fifelse(startsWith(cl, "Class_II"), "DNA",
             fifelse(startsWith(cl, "Class_I"), "oTE",
             fifelse(startsWith(cl, "Simple_repeat") | startsWith(cl, "Low_complexity"), NA_character_, "rep"))))]
  g[!is.na(grp), .(contig, start, end, grp)]
}

rows <- list()
for (sp in sp_all) {
  g <- gff_read(file.path(D, paste0(sp, ".gff"))); te <- te_of(sp)
  cp <- cp_all[species == sp][order(contig)]
  for (i in seq_len(nrow(cp))) {
    r <- cp[i]; flip <- r$strand == "-"
    tox <- function(s, e) if (flip) list(x0 = r$end - e, x1 = r$end - s + 1) else list(x0 = s - r$start, x1 = e - r$start + 1)
    m <- g[type == "mRNA" & contig == r$contig & start == r$start & end == r$end & startsWith(attr1(attr, "Target"), paste0(Q6, " "))][1]
    cds <- g[type == "CDS" & attr1(attr, "Parent") == attr1(m$attr, "ID")]
    ex <- as.data.table(tox(cds$start, cds$end))[order(x0)]
    span <- max(ex$x1)
    intr <- if (nrow(ex) > 1) data.table(a = ex$x1[-nrow(ex)], b = ex$x0[-1]) else data.table(a = numeric(), b = numeric())
    intron_bp <- sum(intr$b - intr$a)
    tt <- NULL; m_te <- list(LTR = NA, oTE = NA, DNA = NA, rep = NA); te_in <- NA; te_up <- NA
    if (!is.null(te)) {
      tt <- te[contig == r$contig]
      if (nrow(tt)) { xx <- tox(tt$start, tt$end); tt[, `:=`(x0 = xx$x0, x1 = xx$x1)]; tt <- tt[x1 > -UP & x0 < span + DOWN] }
      else tt[, `:=`(x0 = numeric(), x1 = numeric())]
      inb <- function(u) sum(mapply(function(a, b) covered(u, a, b), intr$a, intr$b))
      for (k in names(m_te)) { u <- merge_iv(tt[grp == k]$x0, tt[grp == k]$x1); m_te[[k]] <- inb(u) }
      u_te <- merge_iv(tt[grp != "rep"]$x0, tt[grp != "rep"]$x1)
      te_in <- inb(u_te); te_up <- covered(u_te, -UPWIN, 0)
    }
    rows[[length(rows) + 1]] <- list(sp = sp, hap = hap(r$contig), contig = r$contig, strand = r$strand, ex = ex,
      span = span, tt = tt, carp = !is.null(te), intron_bp = intron_bp, te_in = te_in, te_up = te_up, m_te = m_te)
  }
}

# one line per copy
summ <- rbindlist(lapply(rows, function(r) data.table(
  species = r$sp, centromere = if (is.na(cen[r$sp])) "other" else cen[r$sp], contig = r$contig, strand = r$strand,
  n_exons = nrow(r$ex), span_kb = round(r$span / 1000, 2), intron_kb = round(r$intron_bp / 1000, 2),
  carp = r$carp,
  intron_TE_pct = if (r$carp && r$intron_bp > 0) round(100 * r$te_in / r$intron_bp, 1) else NA_real_,
  intron_LTR_kb = round(r$m_te$LTR / 1000, 2), intron_otherTE_kb = round(r$m_te$oTE / 1000, 2),
  intron_DNA_kb = round(r$m_te$DNA / 1000, 2), intron_repeat_notTE_kb = round(r$m_te$rep / 1000, 2),
  upstream_TE_pct = if (r$carp) round(100 * r$te_up / UPWIN, 1) else NA_real_)))
fwrite(summ, file.path(TED, "te_loci_summary.tsv"), sep = "\t")

rna_txt <- function(sp) {
  if (is.null(rna)) return(NULL)
  t6 <- rna[species == sp & set == "TPXL6"]
  if (!nrow(t6)) list(txt = "no RNA", none = TRUE) else list(txt = format(t6$junction_reads[1], big.mark = ","), none = FALSE)
}

n <- length(rows); maxspan <- max(sapply(rows, `[[`, "span"))
X0 <- -UP; X1 <- ceiling((maxspan + DOWN) / 5000) * 5000
u <- (X1 - X0) / TRACK_IN                                     # bp per inch
xc <- X1 + (0.35 + cumsum(COLS_IN)) * u                       # right edge of each column
W_IN <- NAME_IN + TRACK_IN + 0.35 + sum(COLS_IN) + 0.2
H_IN <- 2.25 + 0.36 * n
draw <- function() {
  par(mai = c(0.75, 0.1, 1.6, 0.1), xpd = NA, family = "sans")
  plot(NA, xlim = c(X0 - NAME_IN * u, X1 + (0.35 + sum(COLS_IN) + 0.1) * u), ylim = c(n + 0.6, 0.4), xaxs = "i",
       axes = FALSE, xlab = "", ylab = "")
  at <- c(X0 / 1000, pretty(c(0, X1 / 1000))); at <- unique(at[at * 1000 >= X0 & at * 1000 <= X1])
  axis(1, at = at * 1000, labels = at, cex.axis = 0.72, col = "grey55", col.axis = INK2, line = 0.2)
  mtext("kb from the first coding base (negative = upstream)", side = 1, line = 2.3, at = (X0 + X1) / 2, cex = 0.7, col = INK2)
  segments(0, 0.55, 0, n + 0.45, lty = 3, col = "#b5b4ae")
  hdr <- c("introns", "TE in\nintrons", sprintf("TE in %d kb\nupstream", UPWIN / 1000), "leaf RNA,\nspliced reads")
  text(xc, 0.05, hdr, adj = c(1, 0.5), cex = 0.66, font = 2, col = INK)
  sps <- sapply(rows, `[[`, "sp"); last <- ""
  x_name <- X0 - (NAME_IN - 0.18) * u; x_hap <- X0 - 0.55 * u
  for (k in seq_len(n)) {
    r <- rows[[k]]
    if (r$sp != last) {
      if (k > 1) segments(X0 - (NAME_IN - 0.05) * u, k - 0.5, xc[4], k - 0.5, col = "grey90")
      cc <- if (is.na(cen[r$sp])) NA else CENCOL[cen[r$sp]]
      if (!is.na(cc)) points(X0 - (NAME_IN - 0.08) * u, k, pch = 15, cex = 0.95, col = cc)
      text(x_name, k, nice(r$sp), adj = 0, font = 3, cex = 0.8, col = INK)
      ks <- which(sps == r$sp); yc <- mean(range(ks)); rt <- rna_txt(r$sp)
      if (!is.null(rt)) {
        if (length(ks) > 1) {
          xb <- xc[4] - 0.62 * u; yb <- c(min(ks) - 0.3, max(ks) + 0.3)
          segments(xb, yb[1], xb, yb[2], col = "grey65", lwd = 0.9); segments(xb, yb, xb - 0.06 * u, yb, col = "grey65", lwd = 0.9)
        }
        text(xc[4], yc, rt$txt, adj = c(1, 0.5), cex = if (rt$none) 0.66 else 0.78, font = if (rt$none) 3 else 2,
             col = if (rt$none) MUTED else INK)
      }
      if (!r$carp) {
        why <- st[species == r$sp]$note
        text(xc[3], yc, if (length(why) && startsWith(why[1], "CARP ran on")) "TE annotation is of another assembly" else "no TE annotation yet",
             adj = c(1, 0.5), cex = 0.66, font = 3, col = MUTED)
      }
      last <- r$sp
    }
    text(x_hap, k, r$hap, adj = 0, cex = 0.7, col = INK2)
    segments(-UP, k, r$span + DOWN, k, col = FLANK, lwd = 1)
    segments(0, k, r$span, k, col = LINE, lwd = 1.2)
    if (!is.null(r$tt) && nrow(r$tt)) {
      tt <- r$tt[order(grp != "rep")]                               # grey repeats first, TEs on top
      rect(pmax(tt$x0, -UP), k - 0.21, pmin(tt$x1, r$span + DOWN), k + 0.21, col = GCOL[tt$grp], border = NA)
    }
    rect(r$ex$x0, k - 0.33, pmax(r$ex$x1, r$ex$x0 + MINEX), k + 0.33, col = INK, border = NA)
    text(xc[1], k, sprintf("%.1f kb", r$intron_bp / 1000), adj = c(1, 0.5), cex = 0.72, col = INK)
    if (r$carp) {
      text(xc[2], k, if (r$intron_bp > 0) sprintf("%.0f%%", 100 * r$te_in / r$intron_bp) else "-", adj = c(1, 0.5), cex = 0.72, col = INK)
      text(xc[3], k, sprintf("%.0f%%", 100 * r$te_up / UPWIN), adj = c(1, 0.5), cex = 0.72, col = INK)
    }
  }
  mtext("Transposable elements in and around TPXL6", side = 3, line = 5.45, at = x_name, adj = 0, cex = 1.05, font = 2, col = INK)
  mtext(sprintf("CARP repeat annotation over each TPXL6 copy found by miniprot, 5' to 3', %d kb upstream to %d kb past the last coding base; black = coding exons (drawn at least %d bp wide).",
                UP / 1000, DOWN / 1000, MINEX), side = 3, line = 4.3, at = x_name, adj = 0, cex = 0.68, col = INK2)
  mtext("x = 0 is the first coding base found (exon 2 for a copy whose first exon was not found). TE % counts LTR, other retrotransposons and DNA transposons, not satellites or unclassified repeats.",
        side = 3, line = 3.4, at = x_name, adj = 0, cex = 0.68, col = INK2)
  leg <- function(y, labs, cols) {                           # one legend row: square, label, gap
    x <- X0
    for (i in seq_along(labs)) {
      points(x, y, pch = 15, cex = 1.5, col = cols[i])
      text(x + 0.11 * u, y, labs[i], adj = c(0, 0.5), cex = 0.68, col = INK)
      x <- x + 0.11 * u + strwidth(labs[i], cex = 0.68) + 0.32 * u
    }
  }
  leg(-1.0, GRP[c("LTR", "oTE", "DNA")], GCOL[c("LTR", "oTE", "DNA")])
  leg(-0.42, c(GRP[["rep"]], "coding exon", "holocentric", "monocentric"), c(GCOL[["rep"]], INK, CENCOL))
}
dev_png <- function(f, w, h) if (capabilities("cairo")) png(f, width = w, height = h, units = "in", res = 220, type = "cairo") else png(f, width = w, height = h, units = "in", res = 220)
dev_png(paste0(OUTP, ".png"), W_IN, H_IN); draw(); invisible(dev.off())
pdf(paste0(OUTP, ".pdf"), W_IN, H_IN); draw(); invisible(dev.off())
pr <- copy(summ)[, species := sub("^Drosera_", "D. ", species)]
print(pr[, .(species, contig, intron_kb, carp, intron_TE_pct, intron_LTR_kb, intron_otherTE_kb, intron_DNA_kb,
             intron_repeat_notTE_kb, upstream_TE_pct)], row.names = FALSE)
cat(sprintf("wrote %s.png / .pdf (%.1f x %.1f in) and %s\n", OUTP, W_IN, H_IN, file.path(TED, "te_loci_summary.tsv")))
