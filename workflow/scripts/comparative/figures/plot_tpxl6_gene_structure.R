#!/usr/bin/env Rscript
# plot_tpxl6_gene_structure.R -- TPXL6 gene structure in every genome, drawn from the miniprot alignments
#
# One row per TPXL6 copy found by tpxl6_miniprot.sh and assigned TPXL6 by tpxl6_miniprot_summary.R
# (tpxl6_miniprot_loci.tsv), drawn 5'->3' in the direction of the protein (minus-strand copies flipped).
# Introns are to scale in kb; exons (30-150 bp) are drawn at least MINEX kb wide so they stay visible.
# Genomes with no TPXL6 copy get a row saying what the query found instead.
#   red X   exon carrying an in-frame stop codon (miniprot StopCodon)
#   red fs  exon carrying a frameshift (miniprot Frameshift)
#   text    the call: complete / which protein stretch is not found / stop / frameshift, identity to D. regia
#   columns when the RNA summary of tpxl6_rna_kmer.sh is given: leaf-RNA reads that span a splice junction of
#           TPXL6 (spliced mRNA) and of its TPXL1/5-group relative in the same genome (the library control).
#           tpxl6_rna_kmer.sh searches the k-mers of every TPXL6 copy of a species together, so there is one
#           count per species: it is drawn once, beside a bracket spanning all of that species' copies.
# The page is made as wide as its longest text, so nothing is cut off at the right edge.
# Centromere type: config/species.csv, with D. regia (holocentric) and D. capensis (monocentric) added because
# they are imported genomes that species.csv does not list.
# usage: Rscript plot_tpxl6_gene_structure.R MINIPROT_DIR SPECIES_CSV OUT_PREFIX [RNA_KMER_SUMMARY_TSV]
suppressPackageStartupMessages(library(data.table)); options(width = 200)
a <- commandArgs(TRUE)
if (length(a) < 3) stop("usage: plot_tpxl6_gene_structure.R MINIPROT_DIR SPECIES_CSV OUT_PREFIX [RNA_TSV]")
D <- a[1]; SPCSV <- a[2]; OUTP <- a[3]; RNA <- if (length(a) >= 4) a[4] else ""
Q6 <- "regia_g4881.t1"; MINEX <- 0.3
CEX <- 0.7          # text size of the rows
NAME_IN <- 1.75     # inches left of the gene start: species name and haplotype
GENE_IN <- 4.8      # inches for the longest gene
GAP_IN <- 0.3       # inches between blocks
COL_IN <- 1.0       # inches per RNA column
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
# leaf-RNA columns of a species: spliced TPXL6 reads and spliced reads of the related gene
rna_of <- function(sp) {
  if (is.null(rna)) return(NULL)
  t6 <- rna[species == sp & set == "TPXL6"]; ct <- rna[species == sp & set == "control"]
  if (!nrow(t6)) return(list(t6 = "no RNA", ct = "", none = TRUE))
  list(t6 = format(t6$junction_reads[1], big.mark = ","), none = FALSE,
       ct = if (nrow(ct)) format(ct$junction_reads[1], big.mark = ",") else "")
}
rows <- list()
for (sp in sp_all) {
  cp <- L[species == sp & query == Q6 & assigned == "TPXL6"]
  if (!nrow(cp)) {
    p <- L[species == sp & query == Q6]
    rows[[length(rows) + 1]] <- list(sp = sp, none = TRUE, rna = NULL,
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
    rows[[length(rows) + 1]] <- list(sp = sp, none = FALSE, hap = hap(r$contig), contig = r$contig, span = r$span_kb, ex = ex,
      rna = if (i == 1) rna_of(sp) else NULL,
      txt = if (sp == "Drosera_regia" && r$identity >= 0.999) sprintf("%s; the query itself", call)
            else sprintf("%s; %.0f%% identical to D. regia", call, 100 * r$identity), bad = r$stop > 0 || r$frameshift > 0)
  }
}

dev_png <- function(f, w, h) if (capabilities("cairo")) png(f, width = w, height = h, units = "in", res = 220, type = "cairo") else png(f, width = w, height = h, units = "in", res = 220)
# width of the longest call text, measured on the device that draws it (+10% for other fonts)
CALL_IN <- local({ f <- tempfile(fileext = ".png"); dev_png(f, 4, 4); par(family = "sans")
  w <- max(strwidth(sapply(Filter(function(r) !r$none, rows), `[[`, "txt"), units = "inches", cex = CEX))
  invisible(dev.off()); unlink(f); 1.1 * w })
HAS_RNA <- !is.null(rna)
X_CALL <- GENE_IN + GAP_IN                               # inches from the gene start
X_C1 <- X_CALL + CALL_IN + GAP_IN + COL_IN               # right edge of the TPXL6 column
X_C2 <- X_C1 + COL_IN + 0.15                             # right edge of the related-gene column
W_IN <- NAME_IN + (if (HAS_RNA) X_C2 else X_CALL + CALL_IN) + 0.25
SUB2 <- if (HAS_RNA) sprintf("Right: leaf-RNA reads (%s M per species) spanning a splice junction of TPXL6 (all copies together) and of its closest relative in the same genome (TPXL1/5 group).",
                             paste(unique(range(round(rna$reads_scanned / 1e6))), collapse = "-")) else ""

draw <- function() {
  n <- length(rows); xmax <- max(c(5, sapply(rows, function(r) if (r$none) 0 else r$span))) * 1.02
  u <- xmax / GENE_IN                                    # kb per inch
  par(mai = c(0.7, 0.1, if (HAS_RNA) 1.2 else 1.0, 0.1), xpd = NA, family = "sans")
  plot(NA, xlim = c(-NAME_IN, W_IN - NAME_IN - 0.2) * u, ylim = c(n + 0.6, 0.2), xaxs = "i",
       axes = FALSE, xlab = "", ylab = "")
  axis(1, at = pretty(c(0, xmax)), cex.axis = 0.72, col = "grey50", col.axis = "grey30", line = 0.2)
  mtext("kb from the first coding base", side = 1, line = 2.3, at = xmax / 2, cex = 0.7, col = "grey30")
  x_name <- (-NAME_IN + 0.08) * u; x_hap <- -0.5 * u; x_right <- (W_IN - NAME_IN - 0.25) * u
  if (HAS_RNA) {
    text(X_C1 * u, 0.3, "TPXL6, all copies", adj = c(1, 0.5), cex = 0.68, font = 2, col = "grey20")
    text(X_C2 * u, 0.3, "related gene", adj = c(1, 0.5), cex = 0.68, font = 2, col = "grey20")
    text(X_C2 * u, -0.15, "leaf RNA, spliced reads", adj = c(1, 0.5), cex = 0.68, col = "grey30")
  }
  sps <- sapply(rows, `[[`, "sp")
  last <- ""
  for (k in seq_len(n)) {
    r <- rows[[k]]; col <- colof(r$sp)
    if (r$sp != last) {
      if (k > 1) segments(x_name - 0.05 * u, k - 0.5, x_right, k - 0.5, col = "grey88")
      text(x_name, k, nice(r$sp), adj = 0, font = 3, cex = 0.8, col = col)
      last <- r$sp
    }
    if (!is.null(r$rna)) {                               # one count per species, centred on its copies
      ks <- which(sps == r$sp); yc <- mean(range(ks))
      if (length(ks) > 1) {
        xb <- (X_C1 - 0.55) * u; yb <- c(min(ks) - 0.3, max(ks) + 0.3)
        segments(xb, yb[1], xb, yb[2], col = "grey60", lwd = 0.9)
        segments(xb, yb, xb - 0.07 * u, yb, col = "grey60", lwd = 0.9)
      }
      text(X_C1 * u, yc, r$rna$t6, adj = c(1, 0.5), cex = if (r$rna$none) 0.68 else 0.78,
           font = if (r$rna$none) 3 else 2, col = if (r$rna$none) "grey55" else col)
      text(X_C2 * u, yc, r$rna$ct, adj = c(1, 0.5), cex = 0.72, col = "grey35")
    }
    if (r$none) { text(0, k, r$txt, adj = 0, cex = CEX, col = "grey35"); next }
    text(x_hap, k, r$hap, adj = 0, cex = 0.7, col = "grey40")
    segments(0, k, r$span, k, col = "grey45", lwd = 1.1)
    ex <- r$ex
    for (j in seq_len(nrow(ex))) {
      x1 <- max(ex$x1[j], ex$x0[j] + MINEX)
      rect(ex$x0[j], k - 0.26, x1, k + 0.26, col = col, border = NA)
      if (ex$stopc[j]) text((ex$x0[j] + x1) / 2, k - 0.42, "X", col = "#d7191c", cex = 0.75, font = 2)
      if (ex$fs[j])    text((ex$x0[j] + x1) / 2, k - 0.42, "fs", col = "#d7191c", cex = 0.6, font = 2)
    }
    text(X_CALL * u, k, r$txt, adj = 0, cex = CEX, col = if (r$bad) "#b2182b" else "grey20")
  }
  top <- if (HAS_RNA) 3.4 else 2.5
  mtext("TPXL6 gene structure in every genome", side = 3, line = top, at = x_name, adj = 0, cex = 1.05, font = 2)
  mtext(sprintf("miniprot alignment of the D. regia TPXL6 protein (%d aa) to each genome. One row per copy, 5' to 3', introns to scale; exons drawn at least %.1f kb wide.", ql, MINEX),
        side = 3, line = top - 1.15, at = x_name, adj = 0, cex = 0.7, col = "grey30")
  if (HAS_RNA) mtext(SUB2, side = 3, line = top - 2.05, at = x_name, adj = 0, cex = 0.7, col = "grey30")
  legend(x = 0, y = -0.15, horiz = TRUE, bty = "n", cex = 0.68, x.intersp = 0.5, xjust = 0, yjust = 0.5,
         legend = c("holocentric   ", "monocentric   ", "other plants   ", "X  premature stop   ", "fs  frameshift"),
         pch = c(15, 15, 15, NA, NA), col = c(COLS[["holocentric"]], COLS[["monocentric"]], COLS[["other"]], "#d7191c", "#d7191c"))
}
H_IN <- (if (HAS_RNA) 1.9 else 1.7) + 0.34 * length(rows)
dev_png(paste0(OUTP, ".png"), W_IN, H_IN); draw(); invisible(dev.off())
pdf(paste0(OUTP, ".pdf"), W_IN, H_IN); draw(); invisible(dev.off())
cat(sprintf("rows drawn: %d   page %.1f x %.1f in\n", length(rows), W_IN, H_IN))
for (r in rows) cat(sprintf("  %-20s %-6s %-62s %s\n", r$sp, if (r$none) "" else r$hap, r$txt,
                            if (is.null(r$rna)) "" else if (r$rna$none) "no leaf RNA" else
                              paste0("leaf RNA, spliced reads, all copies: TPXL6 ", r$rna$t6, ", related gene ", r$rna$ct)))
cat(sprintf("wrote %s.png / .pdf\n", OUTP))
