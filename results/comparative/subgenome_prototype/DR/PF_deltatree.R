#!/usr/bin/env Rscript
# ============================================================================
# PF_deltatree -- phylogeny with the per-gene Delta distribution at each tip
#
# Topology from DR/dating/all13R/tree.trees (the corrected rooting):
#     ( (regia, Dionaea), (((scorpioides, paradoxa), binata), capensis) )
# Terminal branch thickness = centromere type. Dionaea defines the Delta axis,
# so its tip carries the two anchors at -1 and +1 rather than a distribution.
#
# IN   DR/out/DR02_gene_delta.csv
# OUT  presentation_figures/PFdelta_tree.{pdf,png}
#      presentation_figures/PFdelta_split.csv
# ============================================================================
suppressMessages(library(data.table))
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }

A_COL <- "#0F6E56"; B_COL <- "#993C1D"; N_COL <- "#6E6D69"
A_FIL <- "#1D9E75"; B_FIL <- "#D85A30"
HOLO <- c("Drosera_regia","Drosera_scorpioides","Drosera_paradoxa")
LW_HOLO <- 4.2; LW_MONO <- 1.4; LW_INT <- 1.4
BW <- 0.10

hr("0. LOAD")
GD <- "DR/out/DR02_gene_delta.csv"
if (!file.exists(GD)) die("%s not found", GD)
G <- fread(GD)[is.finite(delta)]
cat("  copies:", format(nrow(G), big.mark=","), "\n")
S <- G[, .(copies=.N, frac_A=mean(delta<0), median=median(delta),
           strong=mean(abs(delta)>0.5)), by=genome]
S[, ratio := sprintf("%.2f : 1", frac_A/(1-frac_A))]
print(S[order(-frac_A)], row.names=FALSE)
cat(sprintf("\n  pooled: frac_A %.3f (%.2f : 1)\n",
    mean(G$delta<0), mean(G$delta<0)/(1-mean(G$delta<0))))
fwrite(S, file.path(OUT,"PFdelta_split.csv"))

# tip y positions, top to bottom
TIP <- c(Dionaea_muscipula=6, Drosera_regia=5, Drosera_capensis=4,
         Drosera_binata=3, Drosera_paradoxa=2, Drosera_scorpioides=1)
if (!all(setdiff(names(TIP),"Dionaea_muscipula") %in% G$genome))
  die("species in TIP missing from the data")

PX0 <- 54; PX1 <- 101; HALF <- 0.40
dx <- function(d) PX0 + (d + 1) * (PX1 - PX0)/2

draw <- function(){
  par(mar=c(4.0,0.4,3.6,0.4))
  plot(NA, xlim=c(0,108), ylim=c(0.35,6.75), axes=FALSE, xlab="", ylab="")

  seg <- function(x1,y1,x2,y2,lwd) segments(x1,y1,x2,y2,lwd=lwd,lend=1,col="black")
  seg(0,4.3125,2,4.3125,LW_INT); seg(2,3.125,2,5.5,LW_INT)
  seg(2,5.5,16,5.5,LW_INT);      seg(2,3.125,10,3.125,LW_INT)
  seg(16,5,16,6,LW_INT)
  seg(10,2.25,10,4,LW_INT);      seg(10,2.25,16,2.25,LW_INT)
  seg(16,1.5,16,3,LW_INT);       seg(16,1.5,22,1.5,LW_INT)
  seg(22,1,22,2,LW_INT)
  term <- list(list(x=16,y=6,g="Dionaea_muscipula"), list(x=16,y=5,g="Drosera_regia"),
               list(x=10,y=4,g="Drosera_capensis"),  list(x=16,y=3,g="Drosera_binata"),
               list(x=22,y=2,g="Drosera_paradoxa"),  list(x=22,y=1,g="Drosera_scorpioides"))
  for (t in term)
    seg(t$x, t$y, 30, t$y, if (t$g %in% HOLO) LW_HOLO else LW_MONO)

  for (t in term) {
    nm <- if (t$g=="Dionaea_muscipula") "Dionaea" else sub("Drosera_","D. ", t$g)
    text(32, t$y, nm, adj=0, cex=1.15, font=3) }

  for (t in term) {
    y <- t$y
    if (t$g == "Dionaea_muscipula") {
      segments(dx(-1), y, dx(1), y, col="grey78", lwd=1, lty=3)
      points(c(dx(-1), dx(1)), c(y,y), pch=21, cex=2.0, lwd=1.4, col="white",
             bg=c(A_FIL, B_FIL))
      text(dx(-1), y+0.30, "its own A copy", adj=0.5, cex=0.88, col=A_COL)
      text(dx(1),  y+0.30, "its own B copy", adj=0.5, cex=0.88, col=B_COL)
      text(dx(0),  y, "defines the scale", adj=0.5, cex=0.92, col=N_COL)
      next }
    v <- G[genome==t$g, delta]
    d <- density(v, bw=BW, from=-1, to=1, n=512)
    sy <- d$y/max(d$y)*HALF*2
    li <- d$x<=0; ri <- d$x>=0
    polygon(c(dx(d$x[li]), dx(0)), c(y-sy[li], y),
            col=adjustcolor(A_FIL, alpha.f=0.32), border=NA)
    polygon(c(dx(0), dx(d$x[ri])), c(y, y-sy[ri]),
            col=adjustcolor(B_FIL, alpha.f=0.32), border=NA)
    lines(dx(d$x), y-sy, lwd=1.9, col="black")
    segments(dx(0), y, dx(0), y-HALF*2*1.02, lwd=1.0, col="black")
    f <- mean(v<0)
    text(dx(-1)-1.0, y-0.30, sprintf("%.0f%%", 100*f),     adj=1, cex=1.15, font=2, col=A_COL)
    text(dx(1)+1.0,  y-0.30, sprintf("%.0f%%", 100*(1-f)), adj=0, cex=1.15, font=2, col=B_COL) }

  axis(1, at=dx(c(-1,-0.5,0,0.5,1)), labels=c("-1","-0.5","0","+0.5","+1"),
       cex.axis=1.0, lwd=0.6, pos=0.42)
  text(dx(0), 0.06, expression(Delta), cex=1.25)
  text(dx(-1), 0.06, "pure A", adj=0.5, cex=1.0, font=2, col=A_COL)
  text(dx(1),  0.06, "pure B", adj=0.5, cex=1.0, font=2, col=B_COL)

  segments(2, 6.62, 8, 6.62, lwd=LW_HOLO, lend=1); text(10, 6.62, "holocentric", adj=0, cex=0.95)
  segments(30, 6.62, 36, 6.62, lwd=LW_MONO, lend=1); text(38, 6.62, "monocentric", adj=0, cex=0.95)

  mtext("Every species splits the same way", side=3, adj=0, line=1.9, cex=1.25, font=2)
  mtext(sprintf("per-gene assignment for %s Drosera gene copies",
        format(nrow(G), big.mark=",")), side=3, adj=0, line=0.5, cex=0.92, col=N_COL) }

hr("1. RENDER")
W <- 11.5; H <- 6.4
pdf(file.path(OUT,"PFdelta_tree.pdf"), width=W, height=H, useDingbats=FALSE); draw(); dev.off()
png(file.path(OUT,"PFdelta_tree.png"), width=W*160, height=H*160, res=160); draw(); dev.off()
cat("  wrote PFdelta_tree.pdf / .png\n")

hr("DONE")
print(list.files(OUT, pattern="^PFdelta_tree"))
