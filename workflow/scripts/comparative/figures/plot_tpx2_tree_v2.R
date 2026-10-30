#!/usr/bin/env Rscript
# plot_tpx2_tree_v2.R -- the rebuilt TPX2-family tree (tpx2_tree_v2.sh): summary, report figure, full figure
#
# Reads tpx2/v2/tpx2_v2.treefile as inferred (node labels "SH-aLRT/UFBoot") and roots it on Arabidopsis WVD2
# (AT5G28646) for display only; support values move with their splits when rerooting.
#   TPXL6 clade  = MRCA of the six proteins used as anchors since June (D. paradoxa g980, D. regia g4881,
#                  D. roseana h2tg..g133, Arabidopsis AT5G37478, Nepenthes, Dionaea); everything below is computed
#   Printed: clade support and members, Beta / monocentric genes inside it, where D. binata's reading frame sits,
#            the neighbouring group (TPXL6's closest relatives) and what it contains
#   OUT.png/.pdf       report figure: the TPXL6 clade gene by gene; the rest of the family collapsed into grey
#                      groups named by their Arabidopsis members; UFBoot (%) at the key nodes only
#   OUT_full.png/.pdf  the whole family as a circle: every Arabidopsis gene named, TPXL6 clade highlighted,
#                      dots on nodes with SH-aLRT >= 80 and UFBoot >= 95
#   OUT_clade.tsv      TPXL6-clade members with length and alignment occupancy
# Arabidopsis names: Smertenko et al. 2021, J Exp Bot 72:1034, Table 1.
# If tpx2/miniprot/tpxl6_miniprot_loci.tsv exists, the D. binata genome result is added to the footnote.
# usage: Rscript plot_tpx2_tree_v2.R V2_DIR OUT_PREFIX [MINIPROT_LOCI_TSV]
suppressPackageStartupMessages({ library(ape); library(data.table) })
options(width = 200)
a <- commandArgs(TRUE)
if (length(a) < 2) stop("usage: plot_tpx2_tree_v2.R V2_DIR OUT_PREFIX [MINIPROT_LOCI_TSV]")
V2 <- a[1]; OUTP <- a[2]; MP <- if (length(a) >= 3) a[3] else ""
TREE <- file.path(V2, "tpx2_v2.treefile")
for (f in c(TREE, file.path(V2, "tpx2_v2.faa"))) if (!file.exists(f)) stop("missing ", f)
dir.create(dirname(OUTP), showWarnings = FALSE, recursive = TRUE)

ANCH <- c("paradoxa_g980.t1", "regia_g4881.t1", "roseana_h2tg000072l-g133.t1", "Atha_AT5G37478",
          "Nepenthes__BSYO01000005.1_000631.1", "Dionaea__scaffold16_001058.1")
BIN <- "binata_TPXL6_ORF"; WVD2 <- "Atha_AT5G28646"
ATNAME <- c(AT3G04630 = "WDL1", AT5G28646 = "WVD2", AT1G54460 = "WDL2", AT3G23090 = "WDL3", AT2G35880 = "WDL4",
            AT4G32330 = "WDL5", AT2G25480 = "WDL6", AT1G70950 = "WDL7", AT3G01710 = "WDL8", AT3G26050 = "WDL9",
            AT5G44270 = "TPXL7", AT1G03780 = "TPX2", AT3G01015 = "TPXL1", AT5G15510 = "TPXL5", AT5G37478 = "TPXL6",
            AT4G22860 = "TPXL3", AT4G11990 = "TPXL2", AT5G07170 = "TPXL4", AT5G62240 = "TPXL8", AT1G23060 = "MDP40")
COL <- c(holo = "#1b9e77", mono = "#d95f02", Arabidopsis = "#7570b3", Beta = "#e7298a", Dionaea = "#66a61e",
         Nepenthes = "#e6ab02", other = "grey50")
SPN <- c(regia = "D. regia", paradoxa = "D. paradoxa", roseana = "D. roseana", binata = "D. binata", capensis = "D. capensis",
         aliciae = "D. aliciae", tokaiensis = "D. tokaiensis", filiformis = "D. filiformis", scorpioides = "D. scorpioides")
MONO <- c("binata", "capensis", "aliciae", "tokaiensis", "filiformis")
hr  <- function(t) cat("\n==", t, strrep("=", max(3, 92 - nchar(t))), "\n")
lin <- function(l) ifelse(grepl("^paradoxa|^regia|^roseana|^scorpioides", l), "holo",
                   ifelse(grepl("^aliciae|^binata|^capensis|^filiformis|^tokaiensis", l), "mono",
                   ifelse(grepl("^Atha", l), "Arabidopsis", ifelse(grepl("^Beta", l), "Beta",
                   ifelse(grepl("^Dionaea", l), "Dionaea", ifelse(grepl("^Nepenthes", l), "Nepenthes", "other"))))))
agi  <- function(l) ifelse(grepl("AT[1-5CM]G[0-9]{5}", l, ignore.case = TRUE),
                           toupper(sub(".*(AT[1-5CM]G[0-9]{5}).*", "\\1", l, ignore.case = TRUE)), NA)
atn  <- function(l) { g <- agi(l); ifelse(!is.na(g) & g %in% names(ATNAME), ATNAME[g], g) }
fa_len <- function(f) { x <- readLines(f, warn = FALSE); h <- grepl("^>", x)
  setNames(as.integer(tapply(nchar(gsub("[*[:space:]-]", "", x[!h])), cumsum(h)[!h], sum)), sub("^>(\\S+).*", "\\1", x[h])) }
sup  <- function(tr, n) { s <- tr$node.label[n - Ntip(tr)]; if (is.null(s) || is.na(s) || !nzchar(s)) return(c(sh = NA, uf = NA))
                          v <- suppressWarnings(as.numeric(strsplit(s, "/", fixed = TRUE)[[1]])); c(sh = v[1], uf = v[2]) }
fmt  <- function(s) sprintf("SH-aLRT %s / UFBoot %s", ifelse(is.na(s["sh"]), "-", format(s["sh"])), ifelse(is.na(s["uf"]), "-", format(s["uf"])))
par_of  <- function(tr, v) tr$edge[tr$edge[, 2] == v, 1]
kids_of <- function(tr, v) tr$edge[tr$edge[, 1] == v, 2]
tips_of <- function(tr, v) if (v <= Ntip(tr)) tr$tip.label[v] else extract.clade(tr, v)$tip.label
comp <- function(x) { t <- table(factor(lin(x), levels = names(COL))); t <- t[t > 0]; paste(names(t), t, sep = " ", collapse = ", ") }

# ------------------------------------------------------------------ read, root, find the clade
tr <- read.tree(TREE)
if (!WVD2 %in% tr$tip.label) stop("WVD2 (", WVD2, ") not in the tree")
tr <- root(tr, outgroup = WVD2, resolve.root = TRUE, edgelabel = TRUE)
len <- fa_len(file.path(V2, "tpx2_v2.faa"))[tr$tip.label]
occ <- if (file.exists(file.path(V2, "trimmed_occupancy.tsv"))) fread(file.path(V2, "trimmed_occupancy.tsv")) else NULL
iq  <- if (file.exists(file.path(V2, "tpx2_v2.iqtree"))) readLines(file.path(V2, "tpx2_v2.iqtree")) else character()
MODEL <- sub(".*: ", "", grep("Best-fit model", iq, value = TRUE)[1])
SITES <- suppressWarnings(as.integer(sub(".* with ([0-9]+) amino-acid sites.*", "\\1", grep("amino-acid sites", iq, value = TRUE)[1])))
anc <- intersect(ANCH, tr$tip.label)
if (length(anc) < 2) stop("fewer than 2 anchors in the tree")
n6 <- getMRCA(tr, anc); T6 <- tips_of(tr, n6); s6 <- sup(tr, n6)
p6 <- par_of(tr, n6); sp6 <- sup(tr, p6)
SIS <- setdiff(tips_of(tr, p6), T6)

hr("tree")
cat(sprintf("%d proteins | model %s | %s alignment columns | rooted on WVD2 for display\n", Ntip(tr), MODEL, SITES))
hr("TPXL6 clade")
cat("support:", fmt(s6), "\n")
cl <- data.table(tip = T6, lineage = lin(T6), aa = unname(len[T6]))
if (!is.null(occ)) cl <- merge(cl, occ[, .(tip = id, occupancy = fraction)], by = "tip", all.x = TRUE)
print(cl[order(lineage, tip)], row.names = FALSE)
cat(sprintf("\nmembers: %s | anchors inside: %d of %d | Beta inside: %d | monocentric inside: %d (%s)\n",
            comp(T6), length(anc), length(ANCH), sum(lin(T6) == "Beta"), sum(lin(T6) == "mono"),
            paste(T6[lin(T6) == "mono"], collapse = ", ")))
if (BIN %in% tr$tip.label) {
  pb <- par_of(tr, which(tr$tip.label == BIN)); sb <- setdiff(tips_of(tr, pb), BIN)
  cat(sprintf("%s: inside the TPXL6 clade: %s | sister: %s | support for (binata + sister): %s\n",
              BIN, BIN %in% T6, paste(sb, collapse = ", "), fmt(sup(tr, pb))))
}
dro <- T6[lin(T6) %in% c("holo", "mono")]
if (length(dro) >= 2) { md <- getMRCA(tr, dro)
  cat(sprintf("Drosera genes form one group: %s | %s\n", length(tips_of(tr, md)) == length(dro), fmt(sup(tr, md)))) }
dn <- intersect(c("Dionaea__scaffold16_001058.1", "Nepenthes__BSYO01000005.1_000631.1"), T6)
if (length(dn) == 2) { m <- getMRCA(tr, dn)
  cat(sprintf("Dionaea + Nepenthes pair up: %s | %s\n", length(tips_of(tr, m)) == 2, fmt(sup(tr, m)))) }
hr("neighbouring group (TPXL6's closest relatives)")
cat(sprintf("%d genes | split from TPXL6: %s\nmembers: %s\nArabidopsis: %s\nmonocentric species present: %s\n",
            length(SIS), fmt(sp6), comp(SIS), paste(atn(SIS[lin(SIS) == "Arabidopsis"]), collapse = ", "),
            paste(intersect(MONO, sub("_.*", "", SIS)), collapse = ", ")))

# ------------------------------------------------------------------ collapsed report tree
tip_text <- function(l) {
  sp <- sub("_.*", "", l); id <- sub("^[^_]+_+", "", l); L <- len[l]
  ifelse(l == BIN, sprintf("D. binata: reading frame in the genome, not annotated (no leaf RNA)  %d aa", L),
  ifelse(lin(l) == "Arabidopsis", sprintf("Arabidopsis %s (%s)  %d aa", atn(l), agi(l), L),
  ifelse(lin(l) == "Nepenthes", sprintf("Nepenthes gracilis  %d aa", L),
  ifelse(lin(l) == "Dionaea", sprintf("Dionaea muscipula  %d aa", L),
  ifelse(lin(l) == "Beta", sprintf("Beta vulgaris  %d aa", L),
  ifelse(sp %in% names(SPN), sprintf("%s  %s  %d aa", SPN[sp], sub("\\.t[0-9]+$", "", id), L), l))))))
}
root_n <- Ntip(tr) + 1
path <- n6; v <- n6
while (v != root_n) { v <- par_of(tr, v); path <- c(path, v) }
nd <- node.depth.edgelength(tr)
B <- rbindlist(lapply(2:length(path), function(i) rbindlist(lapply(setdiff(kids_of(tr, path[i]), path[i - 1]), function(k) {
  tk <- tips_of(tr, k); ti <- match(tk, tr$tip.label)
  data.table(node = k, rep = tk[1], n = length(tk), depth = max(nd[ti]) - nd[k], stem = nd[k] - nd[path[i]],
             atha = paste(sort(unique(na.omit(atn(tk[lin(tk) == "Arabidopsis"])))), collapse = ", "),
             beta = sum(lin(tk) == "Beta"), mono_sp = length(intersect(MONO, sub("_.*", "", tk))), sister = (i == 2))
})))) 
B[, label := fifelse(n == 1, paste0(tip_text(rep), fifelse(rep == WVD2, "  (display root)", "")),
              sprintf("%s%s: %d genes", fifelse(nzchar(atha), paste0(atha, " group"), "group without an Arabidopsis gene"),
                      fifelse(sister, " (closest relatives)", ""), n))]
B[sister == TRUE & n > 1, label := sprintf("%s; incl. %d sugar beet and %d of 5 monocentric Drosera species", label, beta, mono_sp)]
drop <- unlist(lapply(B$node, function(k) tips_of(tr, k)[-1]))
td <- drop.tip(tr, drop, trim.internal = TRUE, collapse.singles = TRUE)
for (j in seq_len(nrow(B))) { e <- which(td$edge[, 2] == match(B$rep[j], td$tip.label)); td$edge.length[e] <- B$stem[j] + max(B$depth[j], 1e-4) }

mp_note <- ""
if (nzchar(MP) && file.exists(MP)) {
  m <- fread(MP)[species == "Drosera_binata" & assigned == "TPXL6"]
  if (nrow(m)) mp_note <- sprintf(" D. binata genome (miniprot): %d %s, %s.", nrow(m), if (nrow(m) == 1) "copy" else "copies",
                                  paste(sprintf("%s %s", m$contig, m$call), collapse = ", "))
}
draw_report <- function(full_support = FALSE) {
  nt <- Ntip(td); xmax <- max(node.depth.edgelength(td))
  par(mar = c(5.0, 0.6, 3.6, 0.6), xpd = NA)
  plot.phylo(td, show.tip.label = FALSE, x.lim = c(0, xmax * 2.25), y.lim = c(0.4, nt + 0.6), edge.width = 1.2, edge.color = "grey35")
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  m6 <- getMRCA(td, intersect(T6, td$tip.label)); ty <- pp$yy[match(T6, td$tip.label)]
  rect(pp$xx[m6] - xmax * 0.01, min(ty) - 0.45, xmax * 2.2, max(ty) + 0.45, col = adjustcolor("#fdbf6f", 0.28), border = NA)
  for (j in seq_len(nrow(B))) {
    i <- match(B$rep[j], td$tip.label); x1 <- pp$xx[i]; y <- pp$yy[i]; x0 <- x1 - max(B$depth[j], 1e-4)
    if (B$n[j] > 1) polygon(c(x0, x1, x1), c(y, y + 0.34, y - 0.34), col = "grey85", border = "grey55") else
      points(x1, y, pch = 21, cex = 1.0, bg = COL[lin(B$rep[j])], col = COL[lin(B$rep[j])])
    text(x1 + xmax * 0.02, y, B$label[j], adj = 0, cex = 0.72, col = "grey25")
  }
  for (l in intersect(T6, td$tip.label)) {
    i <- match(l, td$tip.label); col <- COL[lin(l)]
    points(pp$xx[i], pp$yy[i], pch = 21, cex = 1.25, bg = if (l == BIN) "white" else col, col = col, lwd = 1.6)
    text(pp$xx[i] + xmax * 0.02, pp$yy[i], tip_text(l), adj = 0, cex = 0.76, col = "grey10", font = if (l == BIN) 3 else 1)
  }
  lab_node <- function(n, txt) text(pp$xx[n], pp$yy[n], txt, adj = c(1.08, -0.45), cex = 0.68, col = "grey20")
  nodes <- if (full_support) (Ntip(td) + 1):(Ntip(td) + td$Nnode) else
    unique(na.omit(c(m6, par_of(td, m6),
                     if (length(dro) >= 2) getMRCA(td, intersect(dro, td$tip.label)),
                     if (BIN %in% td$tip.label) par_of(td, match(BIN, td$tip.label)))))
  for (n in nodes) { s <- sup(td, n); if (!is.na(s["uf"])) lab_node(n, if (full_support) sprintf("%s/%s", s["sh"], s["uf"]) else s["uf"]) }
  add.scale.bar(x = 0, y = 0.2, cex = 0.65, col = "grey30")
  mtext("TPXL6 (MAP20) clade within the TPX2 family", side = 3, line = 2.1, adj = 0, cex = 1.05, font = 2)
  mtext(sprintf("IQ-TREE, %s, %s columns, %d proteins from 12 genomes; rest of the family collapsed (grey). Numbers: %s. Rooted on WVD2 for display.",
                MODEL, SITES, Ntip(tr), if (full_support) "SH-aLRT/UFBoot at every node" else "UFBoot (%) at key nodes"),
        side = 3, line = 0.9, adj = 0, cex = 0.7, col = "grey30")
  legend(x = xmax * 0.45, y = 0.2, horiz = TRUE, pt.cex = 1.2, bty = "n", cex = 0.66, xjust = 0, yjust = 0.5, x.intersp = 0.6, text.width = NA,
         legend = paste0(c("holocentric", "monocentric", "Arabidopsis", "Dionaea", "Nepenthes", "not annotated"), "    "),
         pch = 21, pt.bg = c(COL[c("holo", "mono", "Arabidopsis", "Dionaea", "Nepenthes")], "white"),
         col = c(COL[c("holo", "mono", "Arabidopsis", "Dionaea", "Nepenthes", "mono")]), pt.lwd = 1.4)
  mono_ann <- setdiff(T6[lin(T6) == "mono"], BIN)
  mtext(sprintf("In the TPXL6 clade: %d sugar beet genes, %d annotated monocentric Drosera genes%s.",
                sum(lin(T6) == "Beta"), length(mono_ann), if (length(mono_ann)) paste0(" (", paste(mono_ann, collapse = ", "), ")") else ""),
        side = 1, line = 2.4, adj = 0, cex = 0.66, col = "grey30")
  mtext(sprintf("Open circle: D. binata reading frame, not in its annotation (no leaf RNA support), added from the genome.%s", mp_note),
        side = 1, line = 3.3, adj = 0, cex = 0.66, col = "grey30")
}
draw_full <- function() {
  t2 <- tr; is_at <- lin(t2$tip.label) == "Arabidopsis"
  t2$tip.label <- ifelse(is_at, atn(tr$tip.label), "")
  e6 <- which.edge(tr, T6)
  ec <- rep("grey55", nrow(tr$edge)); ec[e6] <- "black"; ew <- rep(0.6, nrow(tr$edge)); ew[e6] <- 1.8
  par(mar = c(1, 1, 3.2, 1), xpd = NA)
  plot.phylo(t2, type = "fan", show.tip.label = TRUE, cex = 0.55, label.offset = 0.04, edge.color = ec, edge.width = ew,
             tip.color = COL[lin(tr$tip.label)], font = 2)
  short <- !is.na(len) & len < 250
  tiplabels(pch = ifelse(short, 24, 21), bg = ifelse(tr$tip.label == BIN, "white", COL[lin(tr$tip.label)]),
            col = COL[lin(tr$tip.label)], cex = 0.55)
  S <- t(sapply((Ntip(tr) + 1):(Ntip(tr) + tr$Nnode), function(n) sup(tr, n)))
  good <- which(!is.na(S[, "sh"]) & !is.na(S[, "uf"]) & S[, "sh"] >= 80 & S[, "uf"] >= 95)
  nodelabels(node = Ntip(tr) + good, pch = 16, cex = 0.35, col = "grey15", frame = "none")
  mtext("TPX2 family, all proteins (TPXL6 clade in bold black)", side = 3, line = 1.6, adj = 0, cex = 1, font = 2)
  mtext(sprintf("%d proteins, IQ-TREE %s. Dots: SH-aLRT >= 80 and UFBoot >= 95. Arabidopsis genes named; triangles: < 250 aa; open marker: D. binata reading frame.",
                Ntip(tr), MODEL), side = 3, line = 0.5, adj = 0, cex = 0.62, col = "grey30")
  legend("bottomright", bty = "n", cex = 0.62, legend = c("holocentric", "monocentric", "Arabidopsis", "Beta", "Dionaea", "Nepenthes"),
         pch = 21, pt.bg = COL[1:6], col = COL[1:6])
}
dev_png <- function(f, w, h) if (capabilities("cairo")) png(f, width = w, height = h, units = "in", res = 220, type = "cairo") else png(f, width = w, height = h, units = "in", res = 220)
H_IN <- max(4.5, 1.6 + 0.36 * Ntip(td))
dev_png(paste0(OUTP, ".png"), 11, H_IN); draw_report(); invisible(dev.off())
pdf(paste0(OUTP, ".pdf"), 11, H_IN); draw_report(); invisible(dev.off())
dev_png(paste0(OUTP, "_support.png"), 11, H_IN); draw_report(TRUE); invisible(dev.off())
dev_png(paste0(OUTP, "_full.png"), 11, 11); draw_full(); invisible(dev.off())
pdf(paste0(OUTP, "_full.pdf"), 11, 11); draw_full(); invisible(dev.off())
fwrite(cl, paste0(OUTP, "_clade.tsv"), sep = "\t")
hr("files")
cat(sprintf("   %s\n", paste0(OUTP, c(".png", ".pdf", "_support.png", "_full.png", "_full.pdf", "_clade.tsv"))), sep = "")
