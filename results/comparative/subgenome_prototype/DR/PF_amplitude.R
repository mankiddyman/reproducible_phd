#!/usr/bin/env Rscript
# ============================================================================
# PF_amplitude -- compute T_AB and T_spec, predict Delta, compare to noise
#
# T_AB   = divergence between the A and B ancestral genomes, in dS.
#          Measured as dS between Dionaea's two homeologous copies at a locus.
#          The pairing comes from synteny + fractionation, never from Delta.
# T_spec = divergence between a Drosera lineage and Dionaea, in dS.
#          Measured as dS between a Drosera copy and a Dionaea copy. No labels.
#
# Both are LABEL-FREE, so the predicted amplitude is not circular:
#     Delta(A copy) = T_spec/T_AB - 1     Delta(B copy) = 1 - T_spec/T_AB
#     amplitude = |Delta|, separation = 2 x amplitude
#
# Those numbers appear in DR02_label.R and DR03b_dionaea.R as "T_AB ~ 0.57,
# speciation ~ 0.48 from DR00", but NOTHING in DR00 computes them and no script
# contains either value -- they survive only in comments. This recomputes them.
#
# Then: amplitude vs the measured per-gene noise gives the signal-to-noise at
# one gene, and amplitude vs noise/sqrt(n) gives it at block level.
#
# IN   DR/out/pairwise_ks.csv, DR/locus_meta.tsv, fractionation_by_chrpair.csv,
#      DR/out/DR02_gene_delta.csv
# OUT  presentation_figures/PFamplitude.{pdf,png}, PFamplitude_stats.csv
# ============================================================================
suppressMessages(library(data.table))
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }
A_COL <- "#1D9E75"; A_DK <- "#0F6E56"; B_COL <- "#D85A30"; B_DK <- "#993C1D"
NUL <- "#9C9A92"; N_COL <- "#6E6D69"; HI <- "#BA7517"
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")
NOISE_ROBUST <- 0.41; NOISE_SD <- 0.605; BLOCK_N <- 560

hr("0. INPUTS")
K <- fread("DR/out/pairwise_ks.csv")[is.finite(dS) & dS >= 0 & dS < 5 & codons >= 100]
L <- fread("DR/locus_meta.tsv")
P <- fread("fractionation_by_chrpair.csv")
s1 <- P$retained_more == P$chrA
SIDE <- c(setNames(ifelse(s1,"A","B"), P$chrA), setNames(ifelse(s1,"B","A"), P$chrB))
cat(sprintf("  dS rows after filtering: %s\n", format(nrow(K), big.mark=",")))

hr("1. T_AB -- Dionaea's two homeologs")
dio <- L[genome=="Dionaea_muscipula"]
dio[, side := unname(SIDE[chr])]
dio <- dio[!is.na(side)]
KD <- K[sp1=="Dionaea_muscipula" & sp2=="Dionaea_muscipula"]
sd1 <- dio$side[match(KD$seq1, dio$tip)]; sd2 <- dio$side[match(KD$seq2, dio$tip)]
AB <- KD[!is.na(sd1) & !is.na(sd2) & sd1 != sd2]
cat(sprintf("  A-vs-B Dionaea comparisons: %s\n", format(nrow(AB), big.mark=",")))
if (nrow(AB) < 50) die("too few Dionaea A-vs-B pairs")
T_AB <- median(AB$dS)
bt <- quantile(replicate(1000, median(sample(AB$dS, replace=TRUE))), c(.025,.975))
cat(sprintf("  T_AB = %.4f dS   95%% CI [%.4f, %.4f]   (comment said 0.57)\n",
    T_AB, bt[1], bt[2]))
cat(sprintf("  IQR %.3f-%.3f | mean %.3f\n",
    quantile(AB$dS,.25), quantile(AB$dS,.75), mean(AB$dS)))

hr("2. T_spec -- each Drosera vs Dionaea")
KS <- K[(sp1 %in% DROS & sp2=="Dionaea_muscipula") | (sp2 %in% DROS & sp1=="Dionaea_muscipula")]
KS[, dsp := fifelse(sp1 %in% DROS, sp1, sp2)]
KS[, dtip := fifelse(sp1=="Dionaea_muscipula", seq1, seq2)]
KS[, dside := unname(SIDE[dio$chr[match(dtip, dio$tip)]])]
cat(sprintf("  Drosera-Dionaea comparisons: %s | with a Dionaea side: %s\n",
    format(nrow(KS), big.mark=","), format(sum(!is.na(KS$dside)), big.mark=",")))
cat(sprintf("  median dS to Dio_A %.3f | to Dio_B %.3f  (should differ if labels mean anything)\n",
    median(KS$dS[KS$dside=="A"], na.rm=TRUE), median(KS$dS[KS$dside=="B"], na.rm=TRUE)))
cat("  T_spec must use the Drosera copy's OWN side only -- the other side measures T_AB.\n")
GL <- fread("DR/out/DR02_gene_delta.csv")[is.finite(delta)]
lab <- GL[, .(tip, lab = fifelse(delta < 0, "A", "B"))]
KS[, dros_tip := fifelse(sp1 %in% DROS, seq1, seq2)]
KS[, dros_lab := lab$lab[match(dros_tip, lab$tip)]]
KSS <- KS[!is.na(dside) & !is.na(dros_lab) & dside == dros_lab]
cat(sprintf("  same-side comparisons: %s\n", format(nrow(KSS), big.mark=",")))
cat("  NOTE this step uses Delta-derived labels, so T_spec here is NOT label-free.\n")
cat("  It is a consistency check, not an independent prediction.\n")
KS <- KSS
TS <- KS[, .(n=.N, T_spec=median(dS)), by=dsp]
TS[, ratio := T_spec/T_AB]
TS[, amplitude := abs(ratio - 1)]
TS[, sep := 2*amplitude]
TS[, snr_gene  := amplitude/NOISE_ROBUST]
TS[, snr_block := amplitude/(NOISE_ROBUST/sqrt(BLOCK_N))]
setorder(TS, -amplitude)
print(TS[, .(species=sub("Drosera_","D. ",dsp), n, T_spec=round(T_spec,3),
             ratio=round(ratio,3), amplitude=round(amplitude,3),
             snr_gene=round(snr_gene,2), snr_block=round(snr_block,1))],
      row.names=FALSE)
POOL_T <- median(KS$dS); POOL_A <- abs(POOL_T/T_AB - 1)
cat(sprintf("\n  pooled T_spec = %.4f | amplitude = %.3f | separation = %.3f\n",
    POOL_T, POOL_A, 2*POOL_A))
cat(sprintf("  comment claimed T_spec ~ 0.48 -> amplitude 0.16, separation 0.32\n"))

hr("3. SIGNAL AGAINST NOISE")
cat(sprintf("  per-gene noise: robust sd %.3f | plain sd %.3f\n", NOISE_ROBUST, NOISE_SD))
cat(sprintf("  amplitude / noise at ONE gene   : %.2f  -> a gene is a weak vote\n",
    POOL_A/NOISE_ROBUST))
for (n in c(10, 33, 100, BLOCK_N))
  cat(sprintf("  amplitude / (noise/sqrt(%3d))   : %5.1f SE\n", n, POOL_A/(NOISE_ROBUST/sqrt(n))))
cat(sprintf("\n  expected single-gene accuracy at this amplitude: %.0f%%\n",
    100*pnorm(POOL_A/NOISE_ROBUST)))
G <- fread("DR/out/DR02_gene_delta.csv")[is.finite(delta)]
cat(sprintf("  observed |Delta| median: %.3f | observed A-side fraction %.3f\n",
    median(abs(G$delta)), mean(G$delta < 0)))
OBS <- G[, .(m = median(delta)), by=genome]
cat("\n  observed median Delta per species (should be near -amplitude if mostly A):\n")
print(OBS[order(m)], row.names=FALSE)
fwrite(TS, file.path(OUT,"PFamplitude_stats.csv"))

hr("4. FIGURE")
draw <- function(){
  par(mar=c(5.0,5.6,4.4,2.0))
  xr <- c(-1.25, 1.25)
  x <- seq(xr[1], xr[2], length.out=800)
  plot(NA, xlim=xr, ylim=c(0, 1.62), axes=FALSE, xlab="", ylab="")
  dn <- function(m, s) dnorm(x, m, s)/dnorm(0,0,s)
  polygon(c(x,rev(x)), c(dn(-POOL_A, NOISE_ROBUST), rep(0,length(x))),
          col=adjustcolor(A_COL, alpha.f=0.26), border=NA)
  polygon(c(x,rev(x)), c(dn( POOL_A, NOISE_ROBUST), rep(0,length(x))),
          col=adjustcolor(B_COL, alpha.f=0.26), border=NA)
  lines(x, dn(-POOL_A, NOISE_ROBUST), col=A_DK, lwd=2.2)
  lines(x, dn( POOL_A, NOISE_ROBUST), col=B_DK, lwd=2.2)
  sb <- NOISE_ROBUST/sqrt(BLOCK_N)
  lines(x, dn(-POOL_A, sb)*0.62, col=A_DK, lwd=2.6)
  lines(x, dn( POOL_A, sb)*0.62, col=B_DK, lwd=2.6)
  abline(v=0, col="grey45", lwd=1)
  segments(-POOL_A, 1.17, POOL_A, 1.17, col=HI, lwd=2)
  segments(c(-POOL_A,POOL_A), 1.13, c(-POOL_A,POOL_A), 1.21, col=HI, lwd=2)
  text(0, 1.28, sprintf("separation %.2f", 2*POOL_A), adj=0.5, cex=1.0, font=2, col=HI)
  text(-1.18, 1.02, "one gene", adj=0, cex=1.0, font=2, col=N_COL)
  text(-1.18, 0.86, "wide \u2014 A and B overlap almost completely", adj=0, cex=0.84, col=N_COL)
  text(-1.18, 0.62, sprintf("mean of %d genes", BLOCK_N), adj=0, cex=1.0, font=2, col=N_COL)
  text(-1.18, 0.46, "narrow \u2014 no overlap at all", adj=0, cex=0.84, col=N_COL)
  text(-POOL_A, 1.44, "A copies", adj=0.5, cex=1.0, font=2, col=A_DK)
  text( POOL_A, 1.44, "B copies", adj=0.5, cex=1.0, font=2, col=B_DK)
  axis(1, at=seq(-1,1,0.5), cex.axis=0.95, lwd=0.6)
  mtext(expression(Delta), side=1, line=2.8, cex=1.1)
  mtext("Same signal, two scales", side=3, adj=0, line=2.6, cex=1.25, font=2)
  mtext(sprintf("dS between Dionaea's two copies %.2f | dS from a Drosera copy to its own side %.2f | amplitude %.2f | noise %.2f",
        T_AB, POOL_T, POOL_A, NOISE_ROBUST), side=3, adj=0, line=1.2, cex=0.82, col=N_COL)
  mtext(sprintf("a gene separates at %.1f SE; a block of %d at %.0f SE",
        POOL_A/NOISE_ROBUST, BLOCK_N, POOL_A/sb), side=3, adj=0, line=0.0,
        cex=0.88, col=N_COL) }
W <- 9.2; H <- 5.4
pdf(file.path(OUT,"PFamplitude.pdf"), width=W, height=H, useDingbats=FALSE); draw(); dev.off()
png(file.path(OUT,"PFamplitude.png"), width=W*160, height=H*160, res=160); draw(); dev.off()
cat("  wrote PFamplitude.pdf / .png\n")
