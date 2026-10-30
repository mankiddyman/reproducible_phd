#!/usr/bin/env Rscript
# tpxl6_tree_check.R -- what the two TPX2-family trees say about TPXL6 (AT5G37478; called "TPXL6-short" so far)
#
# Reads the trees exactly as inferred (no re-inference) and prints:
#  1. Family tree (tpx2_family.nwk), rooted on Arabidopsis WVD2 (AT5G28646) so that "clade" and "sister" are defined
#     (plot_tpx2_family.R draws the file's own arbitrary root):
#     - the TPXL6 clade = MRCA of the six anchors used in plot_tpx2_family.R: every tip, its length and lineage
#     - support at the clade node, and at the Dionaea + Nepenthes pairing inside it (species tree: Dionaea goes with Drosera)
#     - Beta and monocentric tips inside the clade (the claim needs none)
#     - Arabidopsis members of the sister group (= what TPXL6's closest long relative is)
#     - short (<250 aa) proteins OUTSIDE the clade: terminal branch length and nearest Arabidopsis tip;
#       a monocentric one whose nearest Arabidopsis tip is AT5G37478 would be a misplaced TPXL6 copy
#  2. Short tree (tpx2_short.treefile), rooted on the long paralogs exactly as plot_tpx2.R does:
#     - support for the short lineage, and where the D. binata ORF tip (PSEUDOGENE in its label) sits
# Support values are the trees' own node labels (ultrafast bootstrap in the short tree); rerooting carries them
# as edge labels, so each value stays with its split.
#
# usage: Rscript tpxl6_tree_check.R tpx2_family.nwk tpx2_family.faa tpx2_short.treefile
suppressMessages(library(ape))
a <- commandArgs(TRUE)
if (length(a) < 3) stop("usage: tpxl6_tree_check.R FAMILY_TREE FAMILY_FAA SHORT_TREE")
FAM <- a[1]; FAA <- a[2]; SHORT <- a[3]
for (f in c(FAM, FAA, SHORT)) if (!file.exists(f)) stop("missing input: ", f)

ANCH  <- c("paradoxa_g980.t1", "regia_g4881.t1", "roseana_h2tg000072l-g133.t1", "Atha_AT5G37478",
           "Nepenthes__BSYO01000005.1_000631.1", "Dionaea__scaffold16_001058.1")
WVD2  <- "Atha_AT5G28646"
LONG  <- c("Atha_AT5G15510", "regia_g4900.t1", "capensis_g4300.t1", "paradoxa_g28615.t1", "roseana_h1tg000004l-g260.t1")
ATNAME <- c(Atha_AT5G37478 = "TPXL6 (MAP20 clade)", Atha_AT5G15510 = "TPXL5", Atha_AT5G44270 = "TPXL7",
            Atha_AT5G28646 = "WVD2", Atha_AT1G03780 = "TPX2")
hr  <- function(t) cat("\n==", t, strrep("=", max(3, 90 - nchar(t))), "\n")
lin <- function(l) ifelse(grepl("^paradoxa|^regia|^roseana", l), "holo",
                   ifelse(grepl("^aliciae|^binata|^capensis|^filiformis|^tokaiensis", l), "mono",
                   ifelse(grepl("^Atha", l), "Arabidopsis", ifelse(grepl("^Beta", l), "Beta",
                   ifelse(grepl("^Dionaea", l), "Dionaea", ifelse(grepl("^Nepenthes", l), "Nepenthes", "other"))))))
fa_len <- function(f) {                       # protein lengths from a FASTA, '*' and whitespace not counted
  x <- readLines(f, warn = FALSE); h <- grepl("^>", x)
  id <- sub("^>(\\S+).*", "\\1", x[h]); g <- cumsum(h)
  L <- tapply(nchar(gsub("[*[:space:]]", "", x[!h])), g[!h], sum)
  setNames(as.integer(L[as.character(seq_along(id))]), id)
}
has_sup <- function(tr) !is.null(tr$node.label) && any(nzchar(tr$node.label) & !is.na(tr$node.label))
sup <- function(tr, n) { if (!has_sup(tr)) return("none in file"); s <- tr$node.label[n - Ntip(tr)]
                         if (is.na(s) || !nzchar(s)) "none at this node" else s }
nm  <- function(x) ifelse(x %in% names(ATNAME), paste0(x, " = ", ATNAME[x]), x)

# ---------------------------------------------------------------- 1. family tree
tr <- read.tree(FAM)
cat("family tree:", FAM, "|", Ntip(tr), "tips | support labels in file:", has_sup(tr), "\n")
if (!WVD2 %in% tr$tip.label) stop(WVD2, " (WVD2) not in the family tree; cannot root it")
tr  <- root(tr, outgroup = WVD2, resolve.root = TRUE, edgelabel = TRUE)
len <- fa_len(FAA)[tr$tip.label]
cat("lengths found for", sum(!is.na(len)), "of", Ntip(tr), "tips\n")
anc <- intersect(ANCH, tr$tip.label)
cat("anchors found:", length(anc), "of", length(ANCH), "\n")
if (length(anc) < 2) stop("fewer than 2 anchors in the tree")
n    <- getMRCA(tr, anc)
tips <- extract.clade(tr, n)$tip.label

hr("TPXL6 clade (family tree, rooted on WVD2)")
cat("support at the clade node:", sup(tr, n), "\n\n")
cl <- data.frame(tip = tips, lineage = lin(tips), aa = unname(len[tips]))
print(cl[order(cl$lineage, cl$tip), ], row.names = FALSE)
cat("\ncounts:", paste(names(table(cl$lineage)), table(cl$lineage), sep = "=", collapse = "  "), "\n")
cat("Beta tips in clade:", sum(cl$lineage == "Beta"), "| monocentric tips in clade:", sum(cl$lineage == "mono"),
    "| non-anchor tips:", length(setdiff(tips, anc)), "\n")
cat("note: holo tips should be 5 if this is one locus per genome (regia 1 + paradoxa hap1/hap2 + roseana hap1/hap2)\n")

dn <- intersect(c("Dionaea__scaffold16_001058.1", "Nepenthes__BSYO01000005.1_000631.1"), tips)
if (length(dn) == 2) {
  m <- getMRCA(tr, dn); mt <- extract.clade(tr, m)$tip.label
  if (length(mt) == 2) cat("\nDionaea + Nepenthes are sisters inside the clade, support", sup(tr, m),
                           "(species tree: Dionaea with Drosera)\n")
  else cat("\nDionaea and Nepenthes are not sisters; smallest group holding both has", length(mt), "tips\n")
  dro <- tips[lin(tips) %in% c("holo", "mono")]
  if (length(dro) >= 2) { md <- getMRCA(tr, dro)
    cat("Drosera tips form their own group:", length(extract.clade(tr, md)$tip.label) == length(dro),
        "| support", sup(tr, md), "\n") }
}

p   <- tr$edge[tr$edge[, 2] == n, 1]
sis <- setdiff(extract.clade(tr, p)$tip.label, tips)
hr("sister group of the TPXL6 clade (rooted on WVD2)")
cat("support at the parent node:", sup(tr, p), "| sister size:", length(sis), "tips\n")
cat("lineages:", paste(names(table(lin(sis))), table(lin(sis)), sep = "=", collapse = "  "), "\n")
cat("Arabidopsis members:", paste(nm(grep("^Atha", sis, value = TRUE)), collapse = "; "), "\n")
if (length(sis) > 60) cat("(large sister group: TPXL6's closest relative is not resolved at this level)\n")

hr("short (<250 aa) proteins OUTSIDE the TPXL6 clade")
D    <- cophenetic(tr)
atha <- grep("^Atha", tr$tip.label, value = TRUE)
so   <- tr$tip.label[!is.na(len) & len < 250 & !(tr$tip.label %in% tips)]
if (length(so)) {
  ti <- match(so, tr$tip.label)
  ot <- data.frame(tip = so, lineage = lin(so), aa = unname(len[so]),
                   terminal_branch = round(tr$edge.length[match(ti, tr$edge[, 2])], 3),
                   nearest_Atha = sapply(so, function(s) { d <- D[s, setdiff(atha, s)]; names(d)[which.min(d)] }),
                   stringsAsFactors = FALSE)
  ot$flag <- ifelse(ot$nearest_Atha == "Atha_AT5G37478" & ot$lineage == "mono", "<< nearest is TPXL6: look", "")
  ot$nearest_Atha <- nm(ot$nearest_Atha)
  print(ot[order(-ot$terminal_branch), ], row.names = FALSE)
  cat("\nmedian terminal branch, all tips:", round(median(tr$edge.length[tr$edge[, 2] <= Ntip(tr)]), 3), "\n")
} else cat("none\n")

# ---------------------------------------------------------------- 2. short tree
hr("short tree (rooted on the long paralogs, as in plot_tpx2.R)")
st <- read.tree(SHORT)
rt <- intersect(LONG, st$tip.label)
cat(SHORT, "|", Ntip(st), "tips | support labels in file:", has_sup(st), "| rooting on", length(rt), "long tips\n")
st <- root(st, outgroup = rt, resolve.root = TRUE, edgelabel = TRUE)
sl <- setdiff(st$tip.label, rt)
ns <- getMRCA(st, sl)
cat("support for the short lineage (everything but the long paralogs):", sup(st, ns), "\n")
a2 <- intersect(ANCH, st$tip.label)
if (length(a2) >= 2) { n2 <- getMRCA(st, a2); t2 <- extract.clade(st, n2)$tip.label
  cat("TPXL6 anchors' clade:", length(t2), "tips, support", sup(st, n2), "\n") }
bi <- grep("PSEUDOGENE|^binata", st$tip.label, value = TRUE)
if (length(bi)) for (b in bi) {
  pb <- st$edge[st$edge[, 2] == match(b, st$tip.label), 1]
  sb <- setdiff(extract.clade(st, pb)$tip.label, b)
  cat("\nbinata tip", b, "\n  inside the TPXL6 anchors' clade:", if (exists("t2")) b %in% t2 else NA,
      "\n  its sister:", paste(sb, collapse = ", "), "\n  support for (binata + sister):", sup(st, pb), "\n")
} else cat("no binata tip in the short tree\n")
cat("\nall short-tree tips:\n"); print(data.frame(tip = st$tip.label, lineage = lin(st$tip.label)), row.names = FALSE)
