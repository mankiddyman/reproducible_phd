#!/usr/bin/env Rscript
# ============================================================================
# PF_delta -- what Delta measures: formula, three topologies on one scale
#
#   Delta = [ d(Gene_Drosera, Dio_A) - d(Gene_Drosera, Dio_B) ] / d(Dio_A, Dio_B)
#
# The Drosera copy's own branch is in both numerator terms and cancels; the
# gene's rate multiplies every path and divides out. Bounded to [-1,+1], so
# the ends are Dionaea's own copies: pure A at -1, pure B at +1.
#
# Trees are drawn TO SCALE: horizontal distance is divergence, so the depth at
# which the Drosera copy joins is what sets Delta. Each tree sits directly
# above its own value on the scale.
#
# SHOWN is +/-0.75, an ordinary strong call -- section 0 prints the real
# quantiles so this can be checked against the data.
#
# OUT presentation_figures/PFdelta_intuition.{pdf,png}
# ============================================================================
suppressMessages(library(data.table))
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

A_COL <- "#0F6E56"; B_COL <- "#993C1D"; N_COL <- "#6E6D69"
A_FIL <- "#1D9E75"; B_FIL <- "#D85A30"
SHOW  <- 0.75                  # the illustrative Delta drawn left and right
TSZ <- 1.05; HSZ <- 1.25

hr("0. WHAT THE REAL DELTAS LOOK LIKE")
GD <- "DR/out/DR02_gene_delta.csv"
if (file.exists(GD)) {
  G <- fread(GD)
  G[, chk := (dRA - dRB)/dAB]
  cat(sprintf("  formula check: max |(dRA-dRB)/dAB - delta_raw| = %.3g over %s rows\n",
      max(abs(G$chk - G$delta_raw), na.rm=TRUE), format(nrow(G), big.mark=",")))
  q <- quantile(G$delta, c(.05,.10,.25,.50,.75,.90,.95), na.rm=TRUE)
  cat("  Delta quantiles:\n"); print(round(q, 3))
  cat(sprintf("  sd %.3f | |Delta| > 0.5 : %.1f%% | |Delta| > 0.75 : %.1f%% | |Delta| < 0.16 : %.1f%%\n",
      sd(G$delta, na.rm=TRUE), 100*mean(abs(G$delta) > 0.5, na.rm=TRUE),
      100*mean(abs(G$delta) > 0.75, na.rm=TRUE), 100*mean(abs(G$delta) < 0.16, na.rm=TRUE)))
  cat(sprintf("  illustrative value drawn: +/-%.2f -- that is the %.0fth/%.0fth percentile\n",
      SHOW, 100*mean(G$delta < -SHOW, na.rm=TRUE), 100*mean(G$delta < SHOW, na.rm=TRUE)))
} else cat("  (", GD, " not found -- drawing anyway)\n", sep="")

AX_Y <- 300; AX_L <- 60; AX_R <- 620
d2x <- function(d) AX_L + (d + 1) * (AX_R - AX_L)/2
TY  <- 100

TR_A <- matrix(c(0,53.5,10,53.5, 10,27,10,80, 10,27,50,27, 10,80,100,80,
                 50,9,50,45, 50,9,87,9, 50,45,100,45,
                 87,0,87,18, 87,0,100,0, 87,18,100,18), ncol=4, byrow=TRUE)
TR_0 <- matrix(c(0,53.5,10,53.5, 10,27,10,80, 10,27,30,27, 10,80,100,80,
                 30,9,30,45, 30,9,50,9, 30,45,100,45,
                 50,0,50,18, 50,0,100,0, 50,18,100,18), ncol=4, byrow=TRUE)
TR_B <- matrix(c(0,46.75,10,46.75, 10,13.5,10,80, 10,13.5,50,13.5, 10,80,100,80,
                 50,0,50,27, 50,0,100,0, 50,27,87,27,
                 87,18,87,36, 87,18,100,18, 87,36,100,36), ncol=4, byrow=TRUE)

GENE <- expression(Gene[italic(Drosera)])
tree <- function(xo, M, tips) {
  segments(M[,1]+xo, M[,2]+TY, M[,3]+xo, M[,4]+TY, lwd=1.8, col="black")
  for (t in tips) text(xo+105, t$y+TY+3.5, t$lab, adj=0, cex=TSZ, col=t$col) }

draw_all <- function(){
  par(mar=c(0.3,0.3,0.3,0.3))
  plot(NA, xlim=c(0,680), ylim=c(400,0), axes=FALSE, xlab="", ylab="")
  text(40, 42, "branch length is divergence \u2014 where the copy joins is what Delta reads",
       adj=0, cex=TSZ, col=N_COL)

  pan <- list(list(d=-SHOW, M=TR_A, col=A_COL, tips=list(
                list(y=0,  lab=expression(Dio[A]), col=A_COL),
                list(y=18, lab=GENE,               col="black"),
                list(y=45, lab=expression(Dio[B]), col=B_COL),
                list(y=80, lab="Nep",              col=N_COL))),
              list(d=0, M=TR_0, col=N_COL, tips=list(
                list(y=0,  lab=expression(Dio[A]), col=A_COL),
                list(y=18, lab=expression(Dio[B]), col=B_COL),
                list(y=45, lab=GENE,               col="black"),
                list(y=80, lab="Nep",              col=N_COL))),
              list(d=SHOW, M=TR_B, col=B_COL, tips=list(
                list(y=0,  lab=expression(Dio[A]), col=A_COL),
                list(y=18, lab=expression(Dio[B]), col=B_COL),
                list(y=36, lab=GENE,               col="black"),
                list(y=80, lab="Nep",              col=N_COL))))

  for (p in pan) {
    xo <- d2x(p$d) - 80
    tree(xo, p$M, p$tips)
    text(d2x(p$d), 218, bquote(Delta == .(sprintf("%+.2f", p$d))),
         adj=0.5, cex=HSZ, font=2, col=p$col)
    segments(d2x(p$d), 230, d2x(p$d), AX_Y-10, col="grey78", lwd=0.9, lty=3) }

  rect(AX_L, AX_Y-5, d2x(0), AX_Y+5, col=adjustcolor(A_FIL, alpha.f=0.22), border=NA)
  rect(d2x(0), AX_Y-5, AX_R, AX_Y+5, col=adjustcolor(B_FIL, alpha.f=0.22), border=NA)
  segments(AX_L, AX_Y, AX_R, AX_Y, col="black", lwd=1.2)
  for (v in c(-1,0,1)) segments(d2x(v), AX_Y-9, d2x(v), AX_Y+9, col="black", lwd=1.2)
  for (p in pan) points(d2x(p$d), AX_Y, pch=21, bg=p$col, col="white", cex=2.3, lwd=1.5)

  text(d2x(-1), AX_Y+28, "\u22121", adj=0.5, cex=TSZ)
  text(d2x(0),  AX_Y+28, "0",       adj=0.5, cex=TSZ)
  text(d2x(1),  AX_Y+28, "+1",      adj=0.5, cex=TSZ)
  text(d2x(-1), AX_Y+54, "pure A", adj=0.5, cex=HSZ, font=2, col=A_COL)
  text(d2x(1),  AX_Y+54, "pure B", adj=0.5, cex=HSZ, font=2, col=B_COL)
  text(d2x(-1), AX_Y+76, "Dionaea's own A copy", adj=0.5, cex=TSZ*0.85, col=N_COL)
  text(d2x(1),  AX_Y+76, "Dionaea's own B copy", adj=0.5, cex=TSZ*0.85, col=N_COL)
  text(d2x(0),  AX_Y+54, "Nepenthes", adj=0.5, cex=TSZ*0.85, col=N_COL)
}

draw_formula <- function(){
  par(mar=c(0,0,0,0)); plot(NA, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")
  text(0.5, 0.46, cex=1.30, expression(Delta == frac(
    d(Gene[italic(Drosera)], Dio[A]) - d(Gene[italic(Drosera)], Dio[B]),
    d(Dio[A], Dio[B])))) }

hr("1. RENDER")
W <- 12; H <- 7.6
for (dv in c("pdf","png")) {
  f <- file.path(OUT, paste0("PFdelta_intuition.", dv))
  if (dv=="pdf") pdf(f, width=W, height=H, useDingbats=FALSE)
  else           png(f, width=W*160, height=H*160, res=160)
  layout(matrix(1:2, nrow=2), heights=c(1, 3.8))
  draw_formula(); draw_all()
  dev.off(); cat("  wrote ", f, "\n", sep="") }

hr("DONE")
print(list.files(OUT, pattern="^PFdelta"))
