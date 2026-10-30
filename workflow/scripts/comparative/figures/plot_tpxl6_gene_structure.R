#!/usr/bin/env Rscript
# plot_tpxl6_gene_structure.R -- TPXL6 gene structure in every genome, drawn from the miniprot alignments
#
# One row per TPXL6 copy found by tpxl6_miniprot.sh and assigned TPXL6 by tpxl6_miniprot_summary.R
# (tpxl6_miniprot_loci.tsv), drawn 5'->3' in the direction of the protein (minus-strand copies flipped).
# Introns are to scale in kb; exons (30-150 bp) are drawn at least MINEX kb wide so they stay visible.
# Genomes with no TPXL6 copy get a row saying what the query found instead.
#   red x   exon carrying an in-frame stop codon (miniprot StopCodon)
#   red fs  exon carrying a frameshift (miniprot Frameshift)
#   right   the call: complete / which protein stretch is not found / stop / frameshift, identity to D. regia,
#           and, when tpxl6_rna_check.tsv is given, leaf-RNA reads and splice junctions for that copy
# Centromere type: config/species.csv, with D. regia (holocentric) and D. capensis (monocentric) added because
# they are imported genomes that species.csv does not list.
# usage: Rscript plot_tpxl6_gene_structure.R MINIPROT_DIR SPECIES_CSV OUT_PREFIX [RNA_TSV]
suppressPackageStartupMessages(library(data.table)); options(width = 200)
a <- commandArgs(TRUE)
if (length(a) < 3) stop("usage: plot_tpxl6_gene_structure.R MINIPROT_DIR SPECIES_CSV OUT_PREFIX [RNA_TSV]")
D <- a[1]; SPCSV <- a[2]; OUTP <- a[3]; RNA <- if (length(a) >= 4) a[4] else ""
Q6 <- "regia_g4881.t1"; MINEX <- 0.3
dir.create(dirname(OUTP), showWarnings = FALSE, recursive = TRUE)

ORDER <- c("Nepenthes_gracilis", "Dionaea_muscipula", "Drosera_regia", "Drosera_paradoxa", "Drosera_roseana",
           "Drosera_scorpioides", "Drosera_binata", "Drosera_capensis", "Drosera_aliciae", "Drosera_tokaiensis",
           "Drosera_filiformis")
NICE <- c(Nepenthes_gracilis = "Nepenthes gracilis", Dionaea_muscipula = "Dionaea muscipula")
cen <- c(Drosera_regia = "holocentric", Drosera_capensis = "monocentric")
if (file.exists(SPCSV)) { s <- fread(SPCSV, select = c("species_id", "centromere")); cen[s$species_id] <- s$centromere }
COLS <- c(holocentric = "#1b9e77", monocentric = "#d95f02", other = "grey45")
colof <- function(sp) if (!is.na(cen[sp]) && cen[sp] %in% names(COLS)) COLS[[cen[sp]]] else COLS[["other"]]
nice <- function(sp) if (sp %in% names(NICE)) NICE[[sp]] else sub("^Drosera_", "D. ", sp)
hap  <- function(ctg) { h <- regmatches(ctg, regexpr("hap[0-9]+|h[0-9]+tg", ctg)); if (length(h)) sub("^h([0-9]+)tg$", "hap\\1", h) else "" }
attr1 <- function(x, k) ifelse(grepl(paste0("(^|;)", k, "="), x), sub(paste0(".*(^|;)", k, "=([^;]*).*"), "\\2", x), NA)

ql <- local({ x <- readLines(file.path(D, "queries.faa")); h <- grepl("^>", x)
  setNames(tapply(nchar(x[!h]), cumsum(h)[!h], sum), sub("^>", "", x[h])) })[[Q6]]
L <- fread(file.path(D, "tpxl6_miniprot_loci.tsv"))
gen <- fread(file.path(D, "genomes.tsv"), header = FALSE, col.names = c("species", "genome"))
rna <- if (nzchar(RNA) && file.exists(RNA)) fread(RNA) else NULL
sp_all <- c(intersect(ORDER, gen$species), setdiff(gen$species, ORDER))

# exons of every TPXL6 copy, from the GFF, oriented 5'->3' along the protein
gff_of <- function(sp) {
  x <- readLines(file.path(D, paste0(sp, ".gff"))); x <- x[!startsWith(x, "#") & nzchar(x)]
  if (!length(x)) return(NULL)
  g <- fread(text = x, header = FALSE, sep = "\t", quote = "",
             col.names = c("contig", "src", "type", "start", "end", "score", "strand", "phase", "attr"))
  g[, `:=`(id = attr1(attr, "ID"), parent = attr1(attr, "Parent"), target = attr1(attr, "Target"),
           stopc = !is.na(attr1(attr, "StopCodon")), fs = !is.na(attr1(attr, "Frameshift")))]
  g
}
rows <- list()
for (sp in sp_all) {
  cp <- L[species == sp & query == Q6 & assigned == "TPXL6"]
  if (!nrow(cp)) {
    p <- L[species == sp & query == Q6]
    rows[[length(rows) + 1]] <- list(sp = sp, none = TRUE,
      txt = if (nrow(p)) sprintf("no TPXL6 copy: the query only finds TPXL1/5-group paralogs (best %.0f%% identical)", 100 * max(p$identity))
            else "no TPXL6 copy: no alignment at all")
    next
  }
  g <- gff_of(sp)
  cp <- cp[order(contig)]
  for (i in seq_len(nrow(cp))) {
    r <- cp[i]
    m <- g[type == "mRNA" & contig == r$contig & start == r$start & end == r$end & startsWith(target, paste0(Q6, " "))][1]
    ex <- g[type == "CDS" & parent == m$id][order(start)]
    flip <- r$strand == "-"
    pos0 <- if (flip) r$end else r$start
    ex[, `:=`(x0 = if (flip) pos0 - end else start - pos0, x1 = if (flip) pos0 - start else end - pos0)]
    ex[, `:=`(x0 = x0 / 1000, x1 = (x1 + 1) / 1000)]
    ex <- ex[order(x0)]
    tr <- tstrsplit(ex$target, " "); ts <- as.integer(tr[[2]]); te <- as.integer(tr[[3]])
    cov <- sort(unique(unlist(mapply(seq, ts, te))))
    miss <- setdiff(seq_len(ql), cov)
    miss_txt <- if (length(miss)) { b <- cumsum(c(1, diff(miss) != 1)); paste(tapply(miss, b, function(v) if (length(v) > 1) paste0(min(v), "-", max(v)) else v), collapse = ", ") } else ""
    call <- if (r$stop > 0) "premature stop codon" else if (r$frameshift > 0) "frameshift" else
            if (length(miss)) sprintf("aa %s of %d not found", miss_txt, ql) else "complete, no frameshift or stop"
    rr <- if (!is.null(rna)) rna[species == sp & class == "TPXL6" & contig == r$contig & start == r$start] else NULL
    rtxt <- if (!is.null(rr) && nrow(rr)) sprintf("   leaf RNA: %d reads, %d/%d junctions", rr$reads[1], rr$introns_with_junction[1], rr$introns[1]) else ""
    rows[[length(rows) + 1]] <- list(sp = sp, none = FALSE, hap = hap(r$contig), contig = r$contig, span = r$span_kb, ex = ex,
      txt = if (sp == "Drosera_regia" && r$identity >= 0.999) sprintf("%s; the query itself%s", call, rtxt)
            else sprintf("%s; %.0f%% identical to D. regia%s", call, 100 * r$identity, rtxt), bad = r$stop > 0 || r$frameshift > 0)
  }
}

draw <- function() {
  n <- length(rows); xmax <- max(c(5, sapply(rows, function(r) if (r$none) 0 else r$span))) * 1.02
  par(mar = c(3.4, 0.5, 4.6, 0.5), xpd = NA, family = "sans")
  plot(NA, xlim = c(-0.42 * xmax, 2.05 * xmax), ylim = c(n + 0.6, 0.2), axes = FALSE, xlab = "", ylab = "")
  axis(1, at = pretty(c(0, xmax)), cex.axis = 0.72, col = "grey50", col.axis = "grey30", line = 0.2)
  mtext("kb from the first coding base", side = 1, line = 2.3, at = xmax / 2, cex = 0.7, col = "grey30")
  last <- ""
  for (k in seq_len(n)) {
    r <- rows[[k]]; col <- colof(r$sp)
    if (r$sp != last) {
      if (k > 1) segments(-0.42 * xmax, k - 0.5, 2.05 * xmax, k - 0.5, col = "grey88")
      text(-0.40 * xmax, k, nice(r$sp), adj = 0, font = 3, cex = 0.8, col = col)
      last <- r$sp
    }
    if (r$none) { text(0, k, r$txt, adj = 0, cex = 0.72, col = "grey35"); next }
    text(-0.10 * xmax, k, r$hap, adj = 0, cex = 0.7, col = "grey40")
    segments(0, k, r$span, k, col = "grey45", lwd = 1.1)
    ex <- r$ex
    for (j in seq_len(nrow(ex))) {
      x1 <- max(ex$x1[j], ex$x0[j] + MINEX)
      rect(ex$x0[j], k - 0.26, x1, k + 0.26, col = col, border = NA)
      if (ex$stopc[j]) text((ex$x0[j] + x1) / 2, k - 0.42, "X", col = "#d7191c", cex = 0.75, font = 2)
      if (ex$fs[j])    text((ex$x0[j] + x1) / 2, k - 0.42, "fs", col = "#d7191c", cex = 0.6, font = 2)
    }
    text(1.02 * xmax, k, r$txt, adj = 0, cex = 0.7, col = if (r$bad) "#b2182b" else "grey20")
  }
  mtext("TPXL6 gene structure in every genome", side = 3, line = 2.5, adj = 0, cex = 1.05, font = 2)
  mtext(sprintf("miniprot alignment of the D. regia TPXL6 protein (%d aa) to each genome. One row per copy, 5' to 3', introns to scale; exons drawn at least %.1f kb wide.", ql, MINEX),
        side = 3, line = 1.3, adj = 0, cex = 0.7, col = "grey30")
  legend(x = 0, y = -0.05, horiz = TRUE, bty = "n", cex = 0.68, x.intersp = 0.5, xjust = 0, yjust = 0.5,
         legend = c("holocentric   ", "monocentric   ", "other plants   ", "X  premature stop   ", "fs  frameshift"),
         pch = c(15, 15, 15, NA, NA), col = c(COLS[["holocentric"]], COLS[["monocentric"]], COLS[["other"]], "#d7191c", "#d7191c"))
}
dev_png <- function(f, w, h) if (capabilities("cairo")) png(f, width = w, height = h, units = "in", res = 220, type = "cairo") else png(f, width = w, height = h, units = "in", res = 220)
H_IN <- 1.5 + 0.34 * length(rows)
dev_png(paste0(OUTP, ".png"), 12, H_IN); draw(); invisible(dev.off())
pdf(paste0(OUTP, ".pdf"), 12, H_IN); draw(); invisible(dev.off())
cat("rows drawn:", length(rows), "\n")
for (r in rows) cat(sprintf("  %-20s %-6s %s\n", r$sp, if (r$none) "" else r$hap, r$txt))
cat(sprintf("wrote %s.png / .pdf\n", OUTP))
