#!/usr/bin/env Rscript
# ============================================================================
# PF_whynot23 -- why a true 2:1 is observed as ~56:44
#
# DR02_label.R's model: an A copy centres at -0.16, a B copy at +0.16, each
# with per-gene sd ~0.32. So the separation equals the noise. Mixed 2:1, the
# tails cross zero and get counted on the wrong side -- and because A copies
# are twice as numerous, more A leaks right than B leaks left.
#
#   observed frac_A = (2/3)(1-e) + (1/3)e,  e = P(Z > 0.32/0.32 / 2) = 0.309
#                   = 0.564
#
# OUT presentation_figures/PFdelta_whynot23.{pdf,png}
# ============================================================================
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
A_COL <- "#0F6E56"; B_COL <- "#993C1D"; N_COL <- "#6E6D69"
A_FIL <- "#1D9E75"; B_FIL <- "#D85A30"

MU <- 0.16; SD <- 0.32; WA <- 2/3
e   <- pnorm(0, MU, SD, lower.tail=FALSE)
obs <- WA*(1-e) + (1-WA)*e
cat(sprintf("\n  component means %+.2f / %+.2f | per-gene sd %.2f\n", -MU, MU, SD))
cat(sprintf("  P(a copy lands on the wrong side of 0) = %.3f\n", e))
cat(sprintf("  true 2:1 (%.3f) is observed as %.3f  ->  %.0f : %.0f\n\n",
    WA, obs, 100*obs, 100*(1-obs)))

x  <- seq(-1.6, 1.6, length.out=1200)
ya <- WA*dnorm(x, -MU, SD); yb <- (1-WA)*dnorm(x, MU, SD)

draw <- function(){
  par(mar=c(4.2,0.6,3.4,0.6))
  plot(NA, xlim=c(-1.25,1.25), ylim=c(0, max(ya)*1.30), axes=FALSE, xlab="", ylab="")
  wr <- x >= 0; wl <- x <= 0
  polygon(c(0,x[wr],max(x[wr])), c(0,ya[wr],0), col=adjustcolor(A_FIL, alpha.f=0.55), border=NA)
  polygon(c(min(x[wl]),x[wl],0), c(0,yb[wl],0), col=adjustcolor(B_FIL, alpha.f=0.55), border=NA)
  polygon(c(x,rev(x)), c(ya,rep(0,length(x))), col=adjustcolor(A_FIL, alpha.f=0.16), border=NA)
  polygon(c(x,rev(x)), c(yb,rep(0,length(x))), col=adjustcolor(B_FIL, alpha.f=0.16), border=NA)
  lines(x, ya, lwd=2.2, col=A_COL); lines(x, yb, lwd=2.2, col=B_COL)
  segments(0, 0, 0, max(ya)*1.10, lwd=1.3, col="black")

  text(-MU, max(ya)*1.055, "A copies", adj=0.5, cex=1.2, font=2, col=A_COL)
  text(-MU, max(ya)*0.96,  "two thirds of genes", adj=0.5, cex=0.92, col=A_COL)
  text( 0.62, max(yb)*1.45, "B copies", adj=0.5, cex=1.2, font=2, col=B_COL)
  text( 0.62, max(yb)*1.22, "one third", adj=0.5, cex=0.92, col=B_COL)

  arrows(0.03, max(ya)*0.30, 0.42, max(ya)*0.30, length=0.07, lwd=1.4, col=A_COL)
  text(0.44, max(ya)*0.30, "A genes counted as B", adj=0, cex=0.95, col=A_COL)
  arrows(-0.03, max(ya)*0.13, -0.42, max(ya)*0.13, length=0.07, lwd=1.4, col=B_COL)
  text(-0.44, max(ya)*0.13, "B counted as A", adj=1, cex=0.95, col=B_COL)

  axis(1, at=c(-1,-0.5,0,0.5,1), labels=c("-1","","0","","+1"), cex.axis=1.0, lwd=0.6)
  text(0, -max(ya)*0.14, expression(Delta), cex=1.2, xpd=NA)
  mtext("The groups overlap, so the ratio reads flatter than it is",
        side=3, adj=0, line=1.7, cex=1.15, font=2)
  mtext(sprintf("a true 2 : 1 is observed as %.0f : %.0f \u2014 the gap between A and B is no bigger than the scatter within each",
        100*obs, 100*(1-obs)), side=3, adj=0, line=0.4, cex=0.88, col=N_COL) }

W <- 9; H <- 4.6
pdf(file.path(OUT,"PFdelta_whynot23.pdf"), width=W, height=H, useDingbats=FALSE); draw(); dev.off()
png(file.path(OUT,"PFdelta_whynot23.png"), width=W*160, height=H*160, res=160); draw(); dev.off()
cat("  wrote PFdelta_whynot23.pdf / .png\n")
