#!/usr/bin/env Rscript
# syn4_tree_ab.R -- SYN4 gene tree: one tip per gene, the two Droseraceae paralogs (a, b) labelled
#
# Draws the SYN4 clade of the June kleisin family tree exactly as it was inferred (no re-inference):
#   - rooted on Arabidopsis SYN1, as in syn4_tree_vs_og.R; SYN1-3 themselves are not drawn
#   - one tip per gene: its longest isoform (the same choice the OrthoFinder input makes)
#   - fragments (flagged by syn4_tree_vs_og.R: <50% of the median Droseraceae kleisin length) are kept, drawn open
#   - a = the Droseraceae sub-clade holding D. regia's full-length gene, b = the other one. Each bracket lists the
#     Droseraceae species absent from that sub-clade, or present only as fragments -- computed here, not typed in
#   - tip markers coloured by centromere type; branch lengths as inferred (substitutions per site)
#   - D. scorpioides is not in this tree; its genes in today's SYN4 orthogroup are listed in a note
#
# Inputs : TREE  arabidopsis_map/kleisin_v2.treefile (IQ-TREE; the matching .log gives the support type)
#          TIPS  syn4_check/tips_syn4_clade.tsv, written by syn4_tree_vs_og.R (v2)
#          OGD   today's OrthoFinder Orthogroups/ directory; REPO repository root (D. scorpioides GFF)
# Outputs: OUT.png, OUT.pdf, OUT_tips.tsv (kept tips and their paralog); tables printed to stdout
# Usage  : Rscript syn4_tree_ab.R TREE TIPS OGD REPO OUT
suppressPackageStartupMessages({ library(ape); library(data.table) })
options(width = 200)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("usage: Rscript syn4_tree_ab.R TREE TIPS ORTHOGROUPS_DIR REPO OUT_PREFIX")
TF <- args[1]; TIPSF <- args[2]; OGD <- args[3]; REPO <- args[4]; OUTP <- args[5]
dir.create(dirname(OUTP), showWarnings = FALSE, recursive = TRUE)

DROS <- c("regia", "paradoxa", "roseana", "scorpioides", "binata", "capensis", "aliciae", "tokaiensis", "filiformis")
HOLO <- c("regia", "paradoxa", "roseana", "scorpioides")
COL  <- c(holo = "#1b9e77", mono = "#d95f02", other = "grey62")
SYN1_ROOT <- "(^|[_|])SYN1$|AT5G05490"
hr <- function(t) cat(sprintf("\n== %s\n", t))

# ---------------------------------------------------------------- tree, rooted on SYN1 as syn4_tree_vs_og.R does
tr <- read.tree(TF)
if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
r1 <- grep(SYN1_ROOT, tr$tip.label, ignore.case = TRUE)
if (!length(r1)) stop("no SYN1 tip to root on in ", TF)
rt <- root(tr, outgroup = r1[1], resolve.root = TRUE, edgelabel = TRUE)

# ---------------------------------------------------------------- SYN4-clade tips, one per gene
T <- fread(TIPSF)
need <- c("tip", "group", "core", "gene", "og", "len", "fragment", "in_syn4")
if (length(setdiff(need, names(T)))) stop("tips table lacks columns: ", paste(setdiff(need, names(T)), collapse = " "))
S <- T[in_syn4 == TRUE]
if (!all(S$tip %in% rt$tip.label)) stop("tips table and tree disagree, e.g. ", paste(head(setdiff(S$tip, rt$tip.label), 3), collapse = " "))
syn4_node <- getMRCA(rt, S$tip)
sub <- extract.clade(rt, syn4_node)
if (!setequal(sub$tip.label, S$tip)) stop("the clade spanned by the in_syn4 tips holds other tips too -- tips table is from another tree?")
S[, gkey := fifelse(is.na(gene) | gene == "", core, gene)]
setorderv(S, c("group", "gkey", "len", "tip"), c(1L, 1L, -1L, 1L), na.last = TRUE)
K <- S[, .SD[1], by = .(group, gkey)]
dropped <- S[!tip %in% K$tip]
n_iso <- nrow(dropped[group %in% c("Dionaea", DROS)])          # extra isoforms of Droseraceae genes
dup   <- dropped[!group %in% c("Dionaea", DROS)]$tip           # duplicate entries of one protein (e.g. AT_SYN4 + Atha_SYN4)
drop_txt <- sprintf("%d isoform tip%s dropped%s", n_iso, if (n_iso == 1) "" else "s",
                    if (length(dup)) sprintf(", plus %d duplicate entr%s (%s)", length(dup), if (length(dup) == 1) "y" else "ies",
                                             paste(dup, collapse = ", ")) else "")
pr <- ladderize(keep.tip(sub, K$tip), right = FALSE)
nt <- Ntip(pr)
root_node <- nt + 1L

lab     <- function(nd) if (nd > nt && !is.null(pr$node.label)) pr$node.label[nd - nt] else ""
supp_of <- function(nd) { s <- lab(nd); if (nzchar(s)) sprintf("  (%s)", s) else "" }
tips_of <- function(nd) if (nd <= nt) pr$tip.label[nd] else extract.clade(pr, nd)$tip.label
kids_of <- function(nd) pr$edge[pr$edge[, 1] == nd, 2]

# ---------------------------------------------------------------- the Droseraceae node and its two sub-clades
dr <- K[group %in% c("Dionaea", DROS)]$tip
if (length(dr) < 2 || !is.monophyletic(pr, dr)) stop("the Droseraceae SYN4 tips do not form one clade in this tree")
dr_node <- getMRCA(pr, dr)
dk <- kids_of(dr_node)
if (length(dk) != 2) stop("the Droseraceae node is not a bifurcation")
anchor <- K[group == "regia" & fragment == FALSE]$tip[1]
if (is.na(anchor)) stop("no full-length D. regia gene to anchor paralog a")
a_node <- dk[vapply(dk, function(k) anchor %in% tips_of(k), TRUE)]
b_node <- setdiff(dk, a_node)
K[, paralog := fifelse(tip %in% tips_of(a_node), "a", fifelse(tip %in% tips_of(b_node), "b", "other plants"))]
oth <- K[paralog == "other plants"]$tip
oth_node <- if (length(oth) >= 2 && is.monophyletic(pr, oth)) getMRCA(pr, oth) else NA_integer_

spp  <- intersect(c("Dionaea", DROS), unique(K$group))
nice <- function(g) ifelse(g %in% DROS, paste0("D. ", g), g)
state_of <- function(s, p) {
  x <- K[group == s & paralog == p]
  if (!nrow(x)) "absent" else if (all(x$fragment)) "fragment only" else "present"
}
st <- data.table(species = spp)
st[, a := vapply(species, state_of, "", p = "a")]
st[, b := vapply(species, state_of, "", p = "b")]
note_for <- function(p) {
  ab <- st[get(p) == "absent"]$species
  fr <- st[get(p) == "fragment only"]$species
  out <- c(if (length(ab)) paste0("absent: ", paste(nice(ab), collapse = ", ")),
           if (length(fr)) paste0("fragments only: ", paste(nice(fr), collapse = ", ")))
  if (!length(out)) "full-length in every Droseraceae species drawn" else out
}

# ---------------------------------------------------------------- D. scorpioides note (not in this tree)
cur_og <- names(sort(table(K[group %in% DROS & fragment == FALSE & !is.na(og) & og != ""]$og), decreasing = TRUE))[1]
sc_note <- "D. scorpioides is not in this tree"
ogf <- file.path(OGD, "Orthogroups.tsv")
if (length(cur_og) && !is.na(cur_og) && file.exists(ogf)) {
  hdr <- strsplit(readLines(ogf, n = 1L), "\t", fixed = TRUE)[[1]]
  j <- grep("scorpioides", hdr, ignore.case = TRUE)[1]
  ln <- suppressWarnings(system2("grep", c("-m1", "-P", shQuote(paste0("^", cur_og, "\\t")), shQuote(ogf)), stdout = TRUE))
  if (!is.na(j) && length(ln)) {
    v  <- strsplit(ln[1], "\t", fixed = TRUE)[[1]]
    v  <- c(v, rep("", max(0, length(hdr) - length(v))))
    sg <- trimws(strsplit(v[j], ",", fixed = TRUE)[[1]]); sg <- sg[nzchar(sg)]
    where <- character()
    gff <- file.path(REPO, "results", "Drosera_scorpioides", "annotation", "final", "Drosera_scorpioides.final.gff3")
    if (length(sg) && file.exists(gff)) {
      g <- tryCatch(fread(cmd = sprintf("awk -F'\\t' '$3==\"gene\"' '%s'", gff), sep = "\t", header = FALSE, quote = "",
                          select = c(1, 9), colClasses = "character"), error = function(e) data.table())
      if (nrow(g)) {
        setnames(g, c("seqid", "attr"))
        g[, id := sub("^(.*;)?ID=([^;]+).*$", "\\2", attr)]
        where <- g[id %in% sg, seqid]
      }
    }
    sc_note <- sprintf("D. scorpioides is not in this tree: %d gene%s in the same orthogroup (%s)%s; paralog not yet placed",
                       length(sg), if (length(sg) == 1) "" else "s", cur_og,
                       if (length(where)) paste0(", on ", paste(where, collapse = " and ")) else "")
  }
}

# ---------------------------------------------------------------- support type, from the IQ-TREE command line
logf <- sub("\\.(treefile|contree|nwk|newick|tre)$", ".log", TF)
cmdl <- if (file.exists(logf)) grep("^Command:", readLines(logf, n = 300L, warn = FALSE), value = TRUE)[1] else NA_character_
cl <- paste0(" ", cmdl, " ")
supp <- if (is.na(cmdl)) "node support (type not found in the IQ-TREE log)" else
  if (grepl(" -(B|bb) ", cl) && grepl(" -alrt ", cl)) "SH-aLRT / ultrafast bootstrap" else
  if (grepl(" -(B|bb) ", cl)) "ultrafast bootstrap" else
  if (grepl(" -alrt ", cl)) "SH-aLRT" else
  if (grepl(" -b ", cl)) "standard bootstrap" else "node support"

# ---------------------------------------------------------------- labels, markers, brackets
med_all  <- median(T[group %in% c("Dionaea", DROS) & !is.na(len)]$len)
full_med <- median(K[group %in% c("Dionaea", DROS) & fragment == FALSE]$len, na.rm = TRUE)
row_of <- match(pr$tip.label, K$tip)
labs <- lapply(seq_len(nt), function(i) {
  r <- K[row_of[i]]
  sp <- if (r$group %in% DROS) paste0("D. ", r$group) else r$group
  id <- sub("_(rnasupp|norna)$", "", r$gkey)
  extra <- if (isTRUE(r$fragment)) sprintf(", fragment, %s aa", r$len) else
           if (!is.na(r$len) && r$len < 0.85 * full_med) sprintf(", %s aa", r$len) else ""
  bquote(italic(.(sp)) ~ .(paste0(id, extra)))
})
grp    <- K$group[row_of]
isfrag <- K$fragment[row_of]
mcol <- ifelse(grp %in% HOLO, COL[["holo"]], ifelse(grp %in% c(DROS, "Dionaea"), COL[["mono"]], COL[["other"]]))

br <- list(list(node = a_node, lines = c(paste0("SYN4-a", supp_of(a_node)), note_for("a"))),
           list(node = b_node, lines = c(paste0("SYN4-b", supp_of(b_node)), note_for("b"))))
if (!is.na(oth_node)) br[[3]] <- list(node = oth_node, lines = paste0("SYN4, other plants", supp_of(oth_node)))

notes <- c(
  sprintf("One tip per gene (its longest isoform): %s. Open circles: fragments (<50%% of the median Droseraceae kleisin length, %s aa).",
          drop_txt, format(round(med_all))),
  sprintf("Tree: %s (IQ-TREE, %s), pruned to the SYN4 clade, not re-inferred. Rooted on Arabidopsis SYN1; SYN1-3 not drawn (SYN4 clade support %s).",
          basename(TF), format(file.mtime(TF), "%Y-%m-%d"), if (nzchar(lab(root_node))) lab(root_node) else "n/a"),
  sc_note)

# ---------------------------------------------------------------- layout, measured on the output device type
CEX_TIP <- 0.72; CEX_BR <- 0.74; CEX_NOTE <- 0.64
measure_in <- function(x, cex) {
  f <- tempfile(fileext = ".png")
  if (capabilities("cairo")) png(f, width = 6, height = 6, units = "in", res = 72, type = "cairo") else pdf(NULL)
  on.exit({ dev.off(); unlink(f) })
  plot.new()
  1.08 * max(vapply(x, function(e) strwidth(e, units = "inches", cex = cex), 0))
}
lab_in  <- measure_in(labs, CEX_TIP)
br_in   <- measure_in(as.list(unlist(lapply(br, `[[`, "lines"))), CEX_BR)
note_in <- measure_in(as.list(notes), CEX_NOTE)
MAR <- c(0.3, 0.3, 2.6, 0.3); LINE_IN <- 0.2
OFF_IN <- 0.13; GAP_IN <- 0.2; TICK_IN <- 0.06; TXT_IN <- 0.08; PAD_IN <- 0.1
side_in <- OFF_IN + lab_in + GAP_IN + TXT_IN + br_in + PAD_IN
W_IN <- max(9, side_in + 3.6 + (MAR[2] + MAR[4]) * LINE_IN, note_in + 0.4)
plot_w  <- W_IN - (MAR[2] + MAR[4]) * LINE_IN
tree_in <- plot_w - side_in
xmax <- max(node.depth.edgelength(pr))
u <- xmax / tree_in
NB <- 5
H_IN <- max(6, 0.21 * (nt + NB) + (MAR[1] + MAR[3]) * LINE_IN + 0.3)
key <- unique(c(root_node, dr_node, a_node, b_node, if (!is.na(oth_node)) oth_node))

draw <- function() {
  par(mar = MAR, xpd = NA)
  plot.phylo(pr, show.tip.label = FALSE, x.lim = c(0, plot_w * u), y.lim = c(1 - NB, nt + 0.5),
             edge.width = 1.1, edge.color = "grey25")
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  xx <- pp$xx; yy <- pp$yy
  inn <- (nt + 1L):(nt + pr$Nnode)
  sl <- if (is.null(pr$node.label)) rep("", pr$Nnode) else pr$node.label
  small <- inn[nzchar(sl) & !(inn %in% key)]
  if (length(small)) text(xx[small], yy[small], sl[small - nt], adj = c(1.15, -0.45), cex = 0.48, col = "grey50")
  kk <- key[nzchar(sl[key - nt])]
  if (length(kk)) text(xx[kk], yy[kk], sl[kk - nt], adj = c(1.15, -0.45), cex = 0.64, font = 2, col = "grey10")
  points(xx[1:nt] + 0.05 * u, yy[1:nt], pch = ifelse(isfrag, 1, 21), bg = mcol,
         col = ifelse(isfrag, mcol, "grey20"), cex = 1.05, lwd = ifelse(isfrag, 1.6, 0.7))
  for (i in seq_len(nt)) text(xx[i] + OFF_IN * u, yy[i], labs[[i]], adj = c(0, 0.5), cex = CEX_TIP, col = "grey10")
  xb <- max(xx[1:nt]) + (OFF_IN + lab_in + GAP_IN) * u
  for (b in br) {
    y <- range(yy[match(tips_of(b$node), pr$tip.label)]) + c(-0.3, 0.3)
    segments(xb, y[1], xb, y[2], lwd = 1.3, col = "grey20")
    segments(xb - TICK_IN * u, y, xb, y, lwd = 1.3, col = "grey20")
    n <- length(b$lines)
    text(xb + TXT_IN * u, mean(y) + (seq_len(n) - (n + 1) / 2) * -0.85, b$lines, adj = c(0, 0.5),
         cex = CEX_BR, font = c(2, rep(1, n - 1)), col = "grey10")
  }
  legend(x = 0, y = 0.2, legend = c("holocentric", "monocentric", "other plants", "full-length gene", "fragment"),
         pch = c(21, 21, 21, 21, 1), pt.bg = c(COL[["holo"]], COL[["mono"]], COL[["other"]], "grey85", NA),
         col = "grey20", pt.cex = 1.05, cex = 0.72, horiz = TRUE, bty = "n", xjust = 0, yjust = 1, text.col = "grey10")
  text(0, c(-1.35, -2.1, -2.85), notes, adj = c(0, 0.5), cex = CEX_NOTE, col = "grey25")
  add.scale.bar(x = 0, y = -3.75, cex = 0.62, col = "grey25")
  mtext("SYN4 gene tree: one tip per gene, Droseraceae paralogs a and b", side = 3, line = 1.2, adj = 0, cex = 0.95, font = 2)
  mtext(sprintf("Tip markers coloured by centromere type. Numbers at nodes: %s. Scale bar: substitutions per site.", supp),
        side = 3, line = 0.2, adj = 0, cex = 0.7, col = "grey30")
}

f_png <- paste0(OUTP, ".png"); f_pdf <- paste0(OUTP, ".pdf")
if (capabilities("cairo")) png(f_png, width = W_IN, height = H_IN, units = "in", res = 220, type = "cairo") else
  png(f_png, width = W_IN, height = H_IN, units = "in", res = 220)
draw(); invisible(dev.off())
pdf(f_pdf, width = W_IN, height = H_IN); draw(); invisible(dev.off())
fwrite(K[, .(tip, species = group, gene = gkey, len, fragment, paralog, og)], paste0(OUTP, "_tips.tsv"), sep = "\t")

# ---------------------------------------------------------------- what was drawn, for checking against the figure
hr(sprintf("tree: %s", TF))
cat(sprintf("   SYN4 clade: %d tips in the tree -> %d genes drawn (%s)\n", nrow(S), nrow(K), drop_txt))
cat(sprintf("   support  SYN4 %s | Droseraceae %s | a %s | b %s | other plants %s   [%s]\n", lab(root_node), lab(dr_node),
            lab(a_node), lab(b_node), if (is.na(oth_node)) "-" else lab(oth_node), supp))
if (!is.na(cmdl)) cat("  ", trimws(cmdl), "\n")
if (nrow(dropped)) cat("   tips dropped:", paste(dropped$tip, collapse = " "), "\n")
hr("per species: tree tips -> genes drawn")
iso <- S[group %in% spp, .(tree_tips = .N), by = group][
       K[group %in% spp, .(genes = sum(!fragment), fragments = sum(fragment)), by = group], on = "group"]
iso[, isoform_tips_dropped := tree_tips - genes - fragments]
print(iso[order(match(group, spp))], row.names = FALSE)
hr("paralog a / b per species (full-length genes; fragments counted separately)")
cnt <- K[group %in% spp, .(a_genes = sum(paralog == "a" & !fragment), a_fragments = sum(paralog == "a" & fragment),
                           b_genes = sum(paralog == "b" & !fragment), b_fragments = sum(paralog == "b" & fragment)), by = group]
cnt[, centromere := fifelse(group %in% HOLO, "holo", "mono")]
print(cnt[order(match(group, spp))], row.names = FALSE)
cat("\n   a:", paste(note_for("a"), collapse = "; "), "\n   b:", paste(note_for("b"), collapse = "; "), "\n  ", sc_note, "\n")
hr("files")
cat(sprintf("   %s\n   %s\n   %s_tips.tsv\n", f_png, f_pdf, OUTP))
