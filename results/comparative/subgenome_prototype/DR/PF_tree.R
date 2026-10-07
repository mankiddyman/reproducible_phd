#!/usr/bin/env Rscript
# ============================================================================
# PF_tree -- the 13-tip species tree, drawn from the real newick
#
# Topology: DR/tree/astral/astral_all.tre  (ASTRAL 5.7.8, 1,232 gene trees)
#
# Each node carries a PIE of the three quartet frequencies: q1 = share of gene
# trees agreeing with the drawn branch, q2/q3 = the two alternative
# resolutions. Same convention as a PhyParts pie, read off ASTRAL's own
# annotations. A near-solid pie = little conflict; a three-way split = heavy
# conflict, and near-equal q2/q3 is what incomplete lineage sorting produces.
#
# NOTE: do NOT present the match between q1 and ASTRAL's coalescent branch
# length as independent confirmation -- ASTRAL derives that branch length from
# q1 by inverting the same formula, so the agreement is circular.
# Plotted with ape so the branching comes from the file, not from hand-placed
# coordinates. q1 (share of gene trees agreeing) is annotated at the nodes
# that carry it, parsed from full_A.tre / full_B.tre.
#
# ASTRAL-Pro (DR/tree/pro/astralpro_nep_rooted.tre), which never collapses
# multi-copy tips, recovers the same topology.
#
# OUT presentation_figures/PFtree_AB.{pdf,png}, PFtree_support.csv
# ============================================================================
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
if (!requireNamespace("ape", quietly=TRUE)) stop("ape not installed in this env")
suppressMessages(library(ape))
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }

A_COL <- "#1D9E75"; B_COL <- "#D85A30"; HI <- "#BA7517"; HI_T <- "#854F0B"; N_COL <- "#6E6D69"
Q1C <- "#1D9E75"; Q2C <- "#B4B2A9"; Q3C <- "#E4E2DA"   # pie: agree, alt 1, alt 2
HOLO <- c("regia","scorpioides","paradoxa")
LW_HOLO <- 4.0; LW_MONO <- 1.5

hr("0. TOPOLOGY")
TF <- "DR/tree/astral/astral_all.tre"
if (!file.exists(TF)) die("%s not found", TF)
tr <- read.tree(TF)
if (is.null(tr)) die("could not parse %s", TF)
tr <- root(tr, outgroup="Nepenthes", resolve.root=TRUE)
tr <- ladderize(tr, right=FALSE)
tr$edge.length <- NULL
cat("  tips:", length(tr$tip.label), "\n  ", paste(tr$tip.label, collapse=", "), "\n")
if (length(tr$tip.label) != 13) die("expected 13 tips, got %d", length(tr$tip.label))

hr("1. SUPPORT")
grab <- function(f){ s <- paste(readLines(f, warn=FALSE), collapse="")
  m <- regmatches(s, gregexpr("\\[[^]]*\\]", s))[[1]]
  do.call(rbind, lapply(m, function(b){
    g <- function(k){ r <- regmatches(b, regexpr(paste0(k,"=[0-9.eE-]+"), b))
                      if (length(r)) as.numeric(sub(paste0(k,"="),"",r)) else NA_real_ }
    data.frame(q1=g("q1"), q2=g("q2"), q3=g("q3"), pp1=g("pp1"), EN=g("EN")) })) }
SA <- grab("DR/tree/astral/full_A.tre"); SA$set <- "A"
SB <- grab("DR/tree/astral/full_B.tre"); SB$set <- "B"
SA$node <- c("parscorp","regiaDio","core","binataclade")
SB$node <- c("parscorp","binataclade","regiaDio","core")
S <- rbind(SA, SB)
print(S[,c("set","node","q1","q2","q3","pp1","EN")], row.names=FALSE)
gv <- function(st,nd,cl="q1") S[[cl]][S$set==st & S$node==nd]
for (st in c("A","B")) if (abs(gv(st,"regiaDio") - 0.45) > 0.03)
  die("node mapping wrong: %s regiaDio q1 = %.3f", st, gv(st,"regiaDio"))
write.csv(S, file.path(OUT,"PFtree_support.csv"), row.names=FALSE)

# map each annotated clade to its node number in the plotted tree
mrca_of <- function(tips) getMRCA(tr, tips)
CLADE <- list(
  list(s="A", n="regiaDio",    tips=c("regia_A","Dionaea_A")),
  list(s="A", n="core",        tips=c("capensis_A","binata_A","paradoxa_A","scorpioides_A")),
  list(s="A", n="binataclade", tips=c("binata_A","paradoxa_A","scorpioides_A")),
  list(s="A", n="parscorp",    tips=c("paradoxa_A","scorpioides_A")),
  list(s="B", n="regiaDio",    tips=c("regia_B","Dionaea_B")),
  list(s="B", n="core",        tips=c("capensis_B","binata_B","paradoxa_B","scorpioides_B")),
  list(s="B", n="binataclade", tips=c("binata_B","paradoxa_B","scorpioides_B")),
  list(s="B", n="parscorp",    tips=c("paradoxa_B","scorpioides_B")))
for (k in seq_along(CLADE)) {
  CLADE[[k]]$node <- mrca_of(CLADE[[k]]$tips)
  CLADE[[k]]$q1   <- gv(CLADE[[k]]$s, CLADE[[k]]$n, "q1")
  CLADE[[k]]$q2   <- gv(CLADE[[k]]$s, CLADE[[k]]$n, "q2")
  CLADE[[k]]$q3   <- gv(CLADE[[k]]$s, CLADE[[k]]$n, "q3")
  if (is.na(CLADE[[k]]$node)) die("clade %s/%s not found in the tree",
                                  CLADE[[k]]$s, CLADE[[k]]$n) }


cat("\n  pie slices per node (q1 agree / q2 / q3):\n")
for (k in CLADE) cat(sprintf("    %s %-12s %.3f / %.3f / %.3f   sum %.3f\n",
    k$s, k$n, k$q1, k$q2, k$q3, k$q1+k$q2+k$q3))

pretty_lab <- function(x){
  sp <- sub("_[AB]$", "", x); sg <- sub("^.*_", "", x)
  ifelse(x=="Nepenthes", "Nepenthes",
    ifelse(sp=="Dionaea", paste("Dionaea", sg), paste0("D. ", sp, " ", sg))) }
tipcol <- ifelse(tr$tip.label=="Nepenthes", N_COL,
          ifelse(grepl("^(regia|Dionaea)_", tr$tip.label), HI_T, "black"))
elw <- rep(LW_MONO, nrow(tr$edge))
for (i in seq_len(nrow(tr$edge))) {
  d <- tr$edge[i,2]
  if (d <= length(tr$tip.label)) {
    sp <- sub("_[AB]$", "", tr$tip.label[d])
    if (sp %in% HOLO) elw[i] <- LW_HOLO } }

hr("2. DRAW")
draw <- function(){
  par(mar=c(1.0,0.6,4.0,0.6))
  plot.phylo(tr, show.tip.label=FALSE, edge.width=elw, edge.color="black",
             no.margin=FALSE, plot=FALSE)
  pp0 <- get("last_plot.phylo", envir=.PlotPhyloEnv)
  xmax <- max(pp0$xx)
  plot.phylo(tr, show.tip.label=FALSE, edge.width=elw, edge.color="black",
             no.margin=FALSE, x.lim=c(0, xmax*1.62))
  pp <- get("last_plot.phylo", envir=.PlotPhyloEnv)
  GAP <- xmax*0.02
  cat(sprintf("  tree x range 0..%.2f | label column starts at %.2f\n", xmax, xmax+GAP))
  for (i in seq_along(tr$tip.label))
    text(pp$xx[i] + GAP, pp$yy[i], pretty_lab(tr$tip.label[i]),
         adj=0, cex=1.05, font=3, col=tipcol[i])
  asp <- diff(par("usr")[1:2]) / diff(par("usr")[3:4]) *
         par("pin")[2] / par("pin")[1]
  pie_at <- function(cx, cy, v, r){
    v <- v/sum(v); a0 <- pi/2
    for (j in seq_along(v)) {
      a1 <- a0 - 2*pi*v[j]
      th <- seq(a0, a1, length.out=64)
      polygon(c(cx, cx + r*asp*cos(th)), c(cy, cy + r*sin(th)),
              col=c(Q1C,Q2C,Q3C)[j], border=NA)
      a0 <- a1 }
    th <- seq(0, 2*pi, length.out=128)
    lines(cx + r*asp*cos(th), cy + r*sin(th), col="grey35", lwd=0.6) }
  for (k in CLADE) {
    big <- k$n == "regiaDio"
    pie_at(pp$xx[k$node], pp$yy[k$node], c(k$q1, k$q2, k$q3), if (big) 0.46 else 0.30)
    if (big) text(pp$xx[k$node] - GAP*1.6, pp$yy[k$node] + 0.62, sprintf("%.2f", k$q1),
                  adj=1, cex=0.95, font=2, col=HI_T) }
  ya <- mean(pp$yy[match(c("regia_A","scorpioides_A"), tr$tip.label)])
  yb <- mean(pp$yy[match(c("regia_B","scorpioides_B"), tr$tip.label)])
  text(xmax*1.56, ya, "A", adj=0.5, cex=1.5, font=2, col=A_COL)
  text(xmax*1.56, yb, "B", adj=0.5, cex=1.5, font=2, col=B_COL)
  legend("bottomleft", bty="n", cex=0.88, horiz=TRUE, inset=c(0.02,0.055),
         legend=c("holocentric","monocentric"), lwd=c(LW_HOLO, LW_MONO), seg.len=2)
  legend("bottomleft", bty="n", cex=0.88, horiz=TRUE, inset=c(0.02,-0.02),
         legend=c("gene trees agreeing","alternative 1","alternative 2"),
         fill=c(Q1C,Q2C,Q3C), border="grey35")
  mtext("A and B were analysed separately and give the same tree",
        side=3, adj=0, line=2.1, cex=1.25, font=2)
  mtext("each pie is what the gene trees say at that branch \u2014 every node localPP = 1",
        side=3, adj=0, line=0.7, cex=0.88, col=N_COL) }

W <- 9.5; H <- 7.0
pdf(file.path(OUT,"PFtree_AB.pdf"), width=W, height=H, useDingbats=FALSE); draw(); dev.off()
png(file.path(OUT,"PFtree_AB.png"), width=W*160, height=H*160, res=160); draw(); dev.off()
cat("  wrote PFtree_AB.pdf / .png\n")
