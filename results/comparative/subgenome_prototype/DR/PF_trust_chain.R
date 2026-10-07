#!/usr/bin/env Rscript
# ============================================================================
# PF_trust_chain -- how reliable is one gene's Delta? Three routes, one answer.
#
# ROUTE 1  CROSS-SPECIES, INDEPENDENT OF DELTA
#   The four-point test on Drosera-to-Drosera dS, which was never used to build
#   Delta. Same-label pairing wins in a fraction f of tests. If each copy is
#   mislabelled independently with probability e, the test is right when both
#   copies are right or both are wrong:
#        f = (1-e)^2 + e^2   ->   e = (1 - sqrt(2f - 1)) / 2
#   Per-copy accuracy 1-e then implies, with the measured noise s,
#        amplitude = s * qnorm(1-e)
#   Nothing here uses a block call.
#
# ROUTE 2  AGREEMENT WITH BLOCK CALLS
#   Fraction of genes whose own Delta sign matches their block. Not independent
#   (the block was called from these genes) but a direct read of the same thing.
#
# ROUTE 3  BLOCK MEDIANS
#   The amplitude read straight off DR02_segments.csv. A measurement of the A
#   and B distributions, not a prediction -- the blocks were called using Delta.
#
# NOISE  0.41 (robust sd), measured from Nepenthes genes whose true Delta is 0
#        and from within-block spread. See PF_noise.R.
#
# CAVEAT the (1-e)^2 + e^2 inversion assumes the two copies in a test are
#   independent draws. Under AAB a species has two A copies and one B, so the
#   pairs are not composition-symmetric. Section 1 reports the pair composition
#   so this can be judged.
#
# IN  presentation_figures/PFcross_stats.csv, DR/out/DR02_gene_delta.csv,
#     DR/out/DR02_propagated_blocks.csv, DR/out/DR02_segments.csv
# OUT presentation_figures/PFtrust.{pdf,png}, PFtrust_stats.csv
# ============================================================================
suppressMessages(library(data.table))
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }
A_COL <- "#1D9E75"; A_DK <- "#0F6E56"; B_COL <- "#D85A30"; B_DK <- "#993C1D"
NUL <- "#9C9A92"; N_COL <- "#6E6D69"; HI <- "#BA7517"
NOISE <- 0.41; BLOCK_N <- 560

hr("1. ROUTE 1 -- CROSS-SPECIES FOUR-POINT")
CF <- "presentation_figures/PFcross_stats.csv"
if (!file.exists(CF)) die("%s not found -- run DR/PF_crossspecies.R first", CF)
CS <- fread(CF)
RE <- CS[set=="real"]
f  <- mean(RE$win); n <- nrow(RE)
cat(sprintf("  tests: %s | same-label pairing wins: %.4f\n", format(n, big.mark=","), f))
cat(sprintf("  95%% CI [%.4f, %.4f]\n", binom.test(sum(RE$win), n, 0.5)$conf.int[1],
    binom.test(sum(RE$win), n, 0.5)$conf.int[2]))
if ("set" %in% names(CS) && any(CS$set=="shuffled"))
  cat(sprintf("  control, labels shuffled: %.4f\n", mean(CS[set=="shuffled", win])))
if (f <= 0.5) die("f <= 0.5, the inversion is undefined")
e1  <- (1 - sqrt(2*f - 1))/2
acc1 <- 1 - e1
amp1 <- NOISE * qnorm(acc1)
cat(sprintf("\n  implied per-copy error e = %.3f  ->  accuracy %.3f\n", e1, acc1))
cat(sprintf("  implied amplitude at noise %.2f: %.3f\n", NOISE, amp1))
cat("\n  pair composition (the caveat -- AAB is not symmetric):\n")
print(RE[, .N, by=pair][order(-N)], row.names=FALSE)

hr("2. ROUTE 2 -- AGREEMENT WITH BLOCK CALLS")
G <- fread("DR/out/DR02_gene_delta.csv")[is.finite(delta)]
P <- fread("DR/out/DR02_propagated_blocks.csv")[label %in% c("A","B")]
gk <- if ("gene" %in% names(P)) "gene" else "tip"
M <- merge(G[, .(genome, chr, gene, delta)],
           P[, .(genome, chr, gene=get(gk), label)], by=c("genome","chr","gene"))
acc2 <- mean(ifelse(M$delta < 0, "A", "B") == M$label)
amp2 <- NOISE * qnorm(acc2)
cat(sprintf("  genes: %s | agreement with their block: %.4f\n",
    format(nrow(M), big.mark=","), acc2))
cat(sprintf("  implied amplitude: %.3f   (not independent -- the block used these genes)\n", amp2))

hr("3. ROUTE 3 -- BLOCK MEDIANS")
S <- fread("DR/out/DR02_segments.csv")[label %in% c("A","B")]
amp3 <- median(abs(S$med_delta))
cat(sprintf("  blocks: %d | |median Delta| : median %.3f | IQR %.3f-%.3f\n",
    nrow(S), amp3, quantile(abs(S$med_delta),.25), quantile(abs(S$med_delta),.75)))
cat(sprintf("  A blocks median %.3f | B blocks median %+.3f\n",
    median(S$med_delta[S$label=="A"]), median(S$med_delta[S$label=="B"])))
acc3 <- pnorm(amp3/NOISE)
cat(sprintf("  implied single-gene accuracy: %.3f\n", acc3))

hr("4. THE THREE TOGETHER")
R <- data.table(
  route = c("cross-species four-point", "agreement with block call", "block medians"),
  independent = c("yes", "no", "no"),
  accuracy = c(acc1, acc2, acc3),
  amplitude = c(amp1, amp2, amp3))
print(R[, .(route, independent, accuracy=round(accuracy,3), amplitude=round(amplitude,3))],
      row.names=FALSE)
cat(sprintf("\n  spread across routes: accuracy %.3f-%.3f | amplitude %.3f-%.3f\n",
    min(R$accuracy), max(R$accuracy), min(R$amplitude), max(R$amplitude)))
if (max(R$amplitude)/min(R$amplitude) > 1.5)
  cat("  >>> THE ROUTES DISAGREE BY MORE THAN 1.5x -- do not present these as concordant.\n",
      "      The independent route inverts f = (1-e)^2 + e^2, which assumes the two\n",
      "      copies in a test are independent draws at the same error rate. Under AAB\n",
      "      (two A copies, one B) that is not the case, so e is probably underestimated\n",
      "      and its amplitude overestimated.\n", sep="")
cat(sprintf("  at a block of %d genes, SE = %.4f, so a block sits %.0f SE from zero\n",
    BLOCK_N, NOISE/sqrt(BLOCK_N), median(R$amplitude)/(NOISE/sqrt(BLOCK_N))))
fwrite(R, file.path(OUT,"PFtrust_stats.csv"))

hr("5. FIGURE")
AMP <- median(R$amplitude)
draw <- function(){
  layout(matrix(1:2, nrow=1), widths=c(1.25, 1))

  par(mar=c(5.0,1.2,4.4,1.2))
  x <- seq(-1.6, 1.6, length.out=900)
  ya <- dnorm(x, -AMP, NOISE); yb <- dnorm(x, AMP, NOISE)
  plot(NA, xlim=c(-1.6,1.6), ylim=c(0, max(ya)*1.34), axes=FALSE, xlab="", ylab="")
  wr <- x >= 0; wl <- x <= 0
  polygon(c(x,rev(x)), c(ya,rep(0,length(x))), col=adjustcolor(A_COL,alpha.f=0.18), border=NA)
  polygon(c(x,rev(x)), c(yb,rep(0,length(x))), col=adjustcolor(B_COL,alpha.f=0.18), border=NA)
  polygon(c(0,x[wr],max(x[wr])), c(0,ya[wr],0), col=adjustcolor(A_COL,alpha.f=0.55), border=NA)
  polygon(c(min(x[wl]),x[wl],0), c(0,yb[wl],0), col=adjustcolor(B_COL,alpha.f=0.55), border=NA)
  lines(x, ya, col=A_DK, lwd=2.4); lines(x, yb, col=B_DK, lwd=2.4)
  segments(0,0,0,max(ya)*1.08, lwd=1.2, col="black")
  text(-AMP, max(ya)*1.16, "A copies", adj=0.5, cex=1.1, font=2, col=A_DK)
  text( AMP, max(ya)*1.16, "B copies", adj=0.5, cex=1.1, font=2, col=B_DK)
  arrows(0.05, max(ya)*0.30, 0.55, max(ya)*0.30, length=0.07, lwd=1.5, col=A_DK)
  text(0.58, max(ya)*0.30, sprintf("%.0f%% of A genes\nvote B", 100*(1-AMP/AMP*pnorm(AMP/NOISE))),
       adj=0, cex=0.9, col=A_DK)
  axis(1, at=seq(-1.5,1.5,0.5), cex.axis=0.95, lwd=0.6)
  mtext(expression(Delta * "  for a single gene"), side=1, line=2.9, cex=1.05)
  mtext("One gene is a weighted coin", side=3, adj=0, line=2.6, cex=1.25, font=2)
  mtext(sprintf("amplitude %.2f against noise %.2f -> %.0f%% of genes land on the right side",
        AMP, NOISE, 100*pnorm(AMP/NOISE)), side=3, adj=0, line=1.2, cex=0.88, col=N_COL)
  mtext("shaded = genes that vote for the wrong subgenome",
        side=3, adj=0, line=0.0, cex=0.85, col=N_COL)

  par(mar=c(5.2,5.6,4.4,5.2))
  o <- R[order(amplitude)]
  xl <- c(0, max(o$amplitude)*1.22)
  plot(NA, xlim=xl, ylim=c(0.4, nrow(o)+0.7), axes=FALSE, xlab="", ylab="")
  segments(NOISE, 0.5, NOISE, nrow(o)+0.6, col=HI, lwd=1.6, lty=2)
  text(NOISE, nrow(o)+0.68, sprintf(" noise %.2f", NOISE), adj=0, cex=0.88, col=HI)
  for (i in seq_len(nrow(o))) {
    ind <- o$independent[i]=="yes"
    cl  <- if (ind) A_COL else NUL
    dk  <- if (ind) A_DK  else N_COL
    segments(0, i, o$amplitude[i], i, col=adjustcolor(cl, alpha.f=0.55), lwd=9, lend=1)
    points(o$amplitude[i], i, pch=21, bg=cl, col="white", cex=2.3, lwd=1.5)
    text(o$amplitude[i], i+0.30, sprintf("amplitude %.2f", o$amplitude[i]),
         adj=0.5, cex=0.92, font=2, col=dk)
    text(o$amplitude[i], i-0.30, sprintf("%.0f%% of genes right", 100*o$accuracy[i]),
         adj=0.5, cex=0.88, col=dk)
    lb <- strsplit(o$route[i], " ")[[1]]
    mtext(paste(lb[1:min(2,length(lb))], collapse=" "), side=2, at=i+0.14, las=1,
          line=0.4, cex=0.92, font=if (ind) 2 else 1)
    if (length(lb) > 2)
      mtext(paste(lb[3:length(lb)], collapse=" "), side=2, at=i-0.16, las=1,
            line=0.4, cex=0.92, font=if (ind) 2 else 1)
    mtext(if (ind) "independent\nof \u0394" else "uses \u0394", side=4, at=i, las=1,
          line=0.3, cex=0.80, col=dk) }
  axis(1, cex.axis=0.95, lwd=0.6)
  mtext(expression("amplitude: how far an A or B block mean sits from 0, in " * Delta),
        side=1, line=3.0, cex=0.95)
  mtext("The routes do not agree", side=3, adj=0, line=2.6, cex=1.25, font=2)
  mtext(sprintf("%.2f to %.2f | the independent route is the outlier, and its inversion assumes 1:1 not 2:1",
        min(o$amplitude), max(o$amplitude)), side=3, adj=0, line=1.2,
        cex=0.84, col=N_COL) }
W <- 13; H <- 5.4
pdf(file.path(OUT,"PFtrust.pdf"), width=W, height=H, useDingbats=FALSE); draw(); dev.off()
png(file.path(OUT,"PFtrust.png"), width=W*160, height=H*160, res=160); draw(); dev.off()
cat("  wrote PFtrust.pdf / .png\n")
