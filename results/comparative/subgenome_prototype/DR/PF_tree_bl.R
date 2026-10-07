#!/usr/bin/env Rscript
# ============================================================================
# PF_tree_bl -- the 13-tip tree WITH branch lengths
#
# Source: DR/tree/pro/astralpro_nep_rooted.tre (ASTRAL-Pro, rooted on
# Nepenthes). Its lengths are in COALESCENT UNITS (generations / effective
# population size), NOT time and NOT substitutions. A short branch means few
# generations relative to Ne, so more gene-tree conflict -- which is why the
# regia+Dionaea branch is short and its quartet support is 0.45.
#
# WHAT THIS CANNOT SHOW
#   Divergence times. Coalescent units are not clock-like, and Drosera
#   substitutes ~25% faster than Dionaea (Dionaea-Nepenthes 0.932 dS vs
#   Drosera-Nepenthes 1.168 dS), so even substitution-unit lengths would
#   mislead. Dates need the MCMCtree run.
#
# The previous cladogram (PF_tree.R) deliberately dropped lengths so that no
# timing could be read off it. This one shows them and says what they are.
#
# OUT presentation_figures/PFtree_bl.{pdf,png}
# ============================================================================
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=200)
if (!requireNamespace("ape", quietly=TRUE)) stop("ape not installed in this env")
suppressMessages(library(ape))
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }
A_COL <- "#1D9E75"; B_COL <- "#D85A30"; HI <- "#BA7517"; HI_T <- "#854F0B"; N_COL <- "#6E6D69"
HOLO <- c("regia","scorpioides","paradoxa")
LW_HOLO <- 4.0; LW_MONO <- 1.5

hr("0. TREE")
TF <- "DR/tree/pro/astralpro_nep_rooted.tre"
if (!file.exists(TF)) die("%s not found", TF)
tr <- read.tree(TF)
if (is.null(tr)) die("could not parse %s", TF)
cat("  tips:", length(tr$tip.label), "\n  ", paste(tr$tip.label, collapse=", "), "\n")
if (is.null(tr$edge.length)) die("this tree has no branch lengths")
tr <- ladderize(tr, right=FALSE)
cat(sprintf("  branch lengths: min %.4f | median %.4f | max %.4f\n",
    min(tr$edge.length), median(tr$edge.length), max(tr$edge.length)))
zl <- sum(tr$edge.length <= 1e-8)
if (zl) cat(sprintf("  !! %d branch(es) of length ~0 -- they will be invisible\n", zl))

nm <- function(tips) getMRCA(tr, tips)
KEY <- list(
  list(t=c("regia_A","Dionaea_A"),   lab="regia + Dionaea  (A)"),
  list(t=c("regia_B","Dionaea_B"),   lab="regia + Dionaea  (B)"),
  list(t=c("capensis_A","binata_A","paradoxa_A","scorpioides_A"), lab="core Drosera (A)"),
  list(t=c("capensis_B","binata_B","paradoxa_B","scorpioides_B"), lab="core Drosera (B)"))
cat("\n  length of the branch SUBTENDING each clade, in coalescent units:\n")
for (k in KEY) {
  nd <- nm(k$t)
  if (is.na(nd)) { cat(sprintf("    %-24s not found\n", k$lab)); next }
  e <- which(tr$edge[,2] == nd)
  cat(sprintf("    %-24s %.4f\n", k$lab,
      if (length(e)) tr$edge.length[e] else NA)) }
cat("  ^ short = few generations relative to Ne = heavy gene-tree conflict\n")

pretty_lab <- function(x){
  sp <- sub("_[AB]$", "", x); sg <- sub("^.*_", "", x)
  ifelse(x=="Nepenthes", "Nepenthes",
    ifelse(sp=="Dionaea", paste("Dionaea", sg), paste0("D. ", sp, " ", sg))) }
tipcol <- ifelse(tr$tip.label=="Nepenthes", N_COL,
          ifelse(grepl("^(regia|Dionaea)_", tr$tip.label), HI_T, "black"))
elw <- rep(LW_MONO, nrow(tr$edge))
for (i in seq_len(nrow(tr$edge))) {
  d <- tr$edge[i,2]
  if (d <= length(tr$tip.label) && sub("_[AB]$","",tr$tip.label[d]) %in% HOLO)
    elw[i] <- LW_HOLO }

hr("1. DRAW")
draw <- function(){
  par(mar=c(4.6,0.6,4.0,0.6))
  plot.phylo(tr, show.tip.label=FALSE, edge.width=elw, plot=FALSE)
  pp0 <- get("last_plot.phylo", envir=.PlotPhyloEnv); xmax <- max(pp0$xx)
  plot.phylo(tr, show.tip.label=FALSE, edge.width=elw, edge.color="black",
             x.lim=c(0, xmax*1.58))
  pp <- get("last_plot.phylo", envir=.PlotPhyloEnv); GAP <- xmax*0.018
  for (i in seq_along(tr$tip.label))
    text(pp$xx[i] + GAP, pp$yy[i], pretty_lab(tr$tip.label[i]),
         adj=0, cex=1.05, font=3, col=tipcol[i])
  for (k in KEY) {
    nd <- nm(k$t); if (is.na(nd)) next
    big <- grepl("regia", k$lab)
    points(pp$xx[nd], pp$yy[nd], pch=21, bg=if (big) HI else "white",
           col="black", cex=if (big) 1.8 else 1.1, lwd=0.8) }
  ya <- mean(pp$yy[match(c("regia_A","scorpioides_A"), tr$tip.label)])
  yb <- mean(pp$yy[match(c("regia_B","scorpioides_B"), tr$tip.label)])
  text(xmax*1.52, ya, "A", adj=0.5, cex=1.5, font=2, col=A_COL)
  text(xmax*1.52, yb, "B", adj=0.5, cex=1.5, font=2, col=B_COL)
  add.scale.bar(length=round(xmax/5, 2), cex=0.9, lwd=1.4)
  legend("bottomright", inset=c(0.02,0.02), bty="n", cex=0.9, horiz=TRUE,
         legend=c("holocentric","monocentric"), lwd=c(LW_HOLO, LW_MONO), seg.len=2)
  mtext("Branch lengths in coalescent units", side=3, adj=0, line=2.4, cex=1.25, font=2)
  mtext("short branch = little time relative to population size = heavy gene-tree conflict",
        side=3, adj=0, line=1.1, cex=0.88, col=N_COL)
  mtext("NOT time: these are not clock-like, and Drosera substitutes ~25% faster than Dionaea",
        side=3, adj=0, line=-0.1, cex=0.84, col=HI) }
W <- 9.8; H <- 7.0
pdf(file.path(OUT,"PFtree_bl.pdf"), width=W, height=H, useDingbats=FALSE); draw(); dev.off()
png(file.path(OUT,"PFtree_bl.png"), width=W*160, height=H*160, res=160); draw(); dev.off()
cat("  wrote PFtree_bl.pdf / .png\n")
