#!/usr/bin/env Rscript
# ============================================================================
# PF_noise -- measure Delta's noise, two ways, instead of assuming it
#
# DR02_label.R PREDICTED sd(Delta) ~ 0.32 from first principles: both numerator
# terms are distances near 0.6 with 15-25% relative error, so their difference
# carries about +-0.14, divided by a denominator near 0.525.
#
# ESTIMATE 1  NEPENTHES -- BIAS
#   Nepenthes diverged before the A/B split, so its true Delta is 0 by topology.
#   Its observed centre therefore measures BIAS directly. Its spread is an
#   UPPER BOUND on the noise, not an estimate of it: Nepenthes branches are
#   longer than Drosera branches, so dS there is noisier per comparison.
#
# ESTIMATE 2  WITHIN SYNTENY BLOCK -- MAGNITUDE, PER SPECIES
#   Genes inside one GENESPACE syntenic block share ancestry by construction
#   (the blocks come from synteny, never from Delta). Their Delta values are
#   therefore replicates of one underlying value, so the within-block spread is
#   that species' own noise. This is the estimate that transfers.
#
# WHAT IT BUYS  with noise s and a true separation of 0.32, a block of n genes
#   separates A from B by 0.32/(s/sqrt(n)) standard errors. Printed per species.
#
# IN   DR/out/DR02_gene_delta.csv, DR/locus_meta.tsv, DR02_propagated_blocks.csv
#      ../genespace/results/syntenicBlock_coordinates.csv (if present)
# OUT  presentation_figures/PFnoise.{pdf,png}, PFnoise_stats.csv
# ============================================================================
suppressMessages(library(data.table))
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
GS  <- file.path(dirname(getwd()), "genespace")
hr  <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }
A_COL <- "#1D9E75"; A_DK <- "#0F6E56"; NUL <- "#9C9A92"; N_COL <- "#6E6D69"; HI <- "#BA7517"
PRED <- 0.32; SEP <- 0.32

hr("0. PREDICTED vs OBSERVED SPREAD")
G <- fread("DR/out/DR02_gene_delta.csv")[is.finite(delta)]
cat(sprintf("  Drosera copies with a Delta: %s\n", format(nrow(G), big.mark=",")))
cat(sprintf("  raw sd(Delta) across all copies: %.3f  (includes the +-0.16 true signal)\n",
    sd(G$delta)))
cat(sprintf("  predicted noise sd from DR02's derivation: %.3f\n", PRED))
cat("  the raw sd is an overestimate of noise -- it contains the signal too.\n")

hr("1. NEPENTHES -- IS THE METHOD BIASED?")
LOC <- fread("DR/locus_meta.tsv")
K   <- fread("DR/out/pairwise_ks.csv")[is.finite(dS) & dS >= 0 & dS < 5 & codons >= 100]
PR  <- fread("fractionation_by_chrpair.csv")
s1  <- PR$retained_more == PR$chrA
SIDE <- c(setNames(ifelse(s1,"A","B"), PR$chrA), setNames(ifelse(s1,"B","A"), PR$chrB))
dio <- LOC[genome=="Dionaea_muscipula"]
dio[, side := unname(SIDE[chr])]
dio <- dio[!is.na(side)]
ax <- dcast(dio[, .(locus, side, tip)], locus ~ side, value.var="tip",
            fun.aggregate=function(z) if (length(z)==1) z else NA_character_)
ax <- ax[!is.na(A) & !is.na(B)]
setnames(ax, c("A","B"), c("DA","DB"))
cat(sprintf("  loci with exactly one Dionaea copy per side: %s\n",
    format(nrow(ax), big.mark=",")))
DD <- rbind(K[, .(locus=anchor, a=seq1, b=seq2, d=dS)],
            K[, .(locus=anchor, a=seq2, b=seq1, d=dS)])
dm <- setNames(DD$d, paste(DD$locus, DD$a, DD$b))
gd <- function(l,x,y) unname(dm[paste(l,x,y)])
nep <- LOC[genome=="Nepenthes_gracilis", .(locus, tip)]
NE  <- merge(nep, ax, by="locus")
NE[, `:=`(dRA = gd(locus, tip, DA), dRB = gd(locus, tip, DB), dAB = gd(locus, DA, DB))]
NE <- NE[is.finite(dRA) & is.finite(dRB) & is.finite(dAB) & dAB >= 0.20]
NE[, delta := (dRA - dRB)/dAB]
if (!nrow(NE)) die("no Nepenthes deltas could be computed")
cat(sprintf("  Nepenthes genes with a computable Delta: %s\n", format(nrow(NE), big.mark=",")))
cat(sprintf("  BIAS  : median %+.4f | mean %+.4f | Wilcoxon p = %.3g\n",
    median(NE$delta), mean(NE$delta),
    suppressWarnings(wilcox.test(NE$delta)$p.value)))
cat(sprintf("  SPREAD: sd %.3f | MAD-based sd %.3f | IQR %.3f\n",
    sd(NE$delta), mad(NE$delta), IQR(NE$delta)))
cat("  ^ an UPPER BOUND on Drosera's noise -- Nepenthes branches are longer.\n")

hr("2. WITHIN SYNTENY BLOCK -- THE TRANSFERABLE ESTIMATE")
BF <- file.path(GS,"results","syntenicBlock_coordinates.csv")
P  <- fread("DR/out/DR02_propagated_blocks.csv")[label %in% c("A","B")]
gk <- if ("gene" %in% names(P)) "gene" else "tip"
M  <- merge(G[, .(genome, chr, gene, mid, delta)],
            P[, .(genome, chr, gene=get(gk), region, label)], by=c("genome","chr","gene"))
M[, mb := mid/1e6]
if (file.exists(BF)) {
  B <- fread(BF)
  cat("  syntenicBlock_coordinates.csv columns:", paste(names(B), collapse=", "), "\n")
} else cat("  (block coordinates not found -- using region x chromosome as the unit)\n")
M[, unit := paste(genome, chr, region)]
W <- M[, .(n=.N, s=sd(delta)), by=.(genome, unit)][n >= 8 & is.finite(s)]
cat(sprintf("  units with >= 8 measured genes: %d\n", nrow(W)))
SP <- W[, .(units=.N, genes=sum(n),
            noise=sqrt(weighted.mean(s^2, n-1))), by=genome]
SP[, se_at_100 := noise/10]
SP[, z_at_100  := SEP/(noise/10)]
SP[, n_for_3SE := ceiling((3*noise/SEP)^2)]
setorder(SP, noise)
print(SP[, .(genome, units, genes, noise=round(noise,3),
             z_at_100=round(z_at_100,1), n_for_3SE)], row.names=FALSE)
POOL <- sqrt(weighted.mean(W$s^2, W$n-1))
cat(sprintf("\n  POOLED within-block noise: %.3f  (predicted %.3f)\n", POOL, PRED))
cat(sprintf("  Nepenthes upper bound:     %.3f\n", mad(NE$delta)))
cat(sprintf("  a block of 100 genes separates A from B by %.1f standard errors\n",
    SEP/(POOL/10)))
cat(sprintf("  genes needed for a 3-SE call: %d\n", ceiling((3*POOL/SEP)^2)))
fwrite(SP, file.path(OUT,"PFnoise_stats.csv"))

hr("3. PER-UNIT SPREAD -- the distribution, not one number")
W[, sp := sub("Drosera_","D. ", genome)]
print(W[, .(units=.N, genes=sum(n),
            sd_min=round(min(s),3), sd_q25=round(quantile(s,.25),3),
            sd_med=round(median(s),3), sd_q75=round(quantile(s,.75),3),
            sd_max=round(max(s),3)), by=sp][order(sd_med)], row.names=FALSE)
cat(sprintf("\n  across all %d units: sd ranges %.2f to %.2f, median %.2f\n",
    nrow(W), min(W$s), max(W$s), median(W$s)))
cat(sprintf("  units with sd below the predicted %.2f: %d (%.0f%%)\n",
    PRED, sum(W$s < PRED), 100*mean(W$s < PRED)))
cat(sprintf("  unit size: median %d genes | range %d-%d\n",
    as.integer(median(W$n)), min(W$n), max(W$n)))
cc <- suppressWarnings(cor(W$n, W$s, method="spearman"))
cat(sprintf("  Spearman(unit size, sd) = %+.3f %s\n", cc,
    if (!is.na(cc) && cc > 0.2) "-- bigger units are noisier, so they span real boundaries" else ""))

hr("4. FIGURE")
NEsd  <- sd(NE$delta); NEmad <- mad(NE$delta)
draw <- function(){
  layout(matrix(1:2, nrow=1), widths=c(1, 1.15))

  par(mar=c(4.8,4.4,4.2,1.4))
  d <- density(NE$delta, from=-1.5, to=1.5, bw=0.10)
  plot(NA, xlim=c(-1.5,1.5), ylim=c(0, max(d$y)*1.30), axes=FALSE, xlab="", ylab="")
  polygon(c(d$x, rev(d$x)), c(d$y, rep(0,length(d$y))),
          col=adjustcolor(NUL, alpha.f=0.40), border=NA)
  lines(d$x, d$y, lwd=2.2, col=N_COL)
  abline(v=0, lwd=1.3, col="black")
  abline(v=median(NE$delta), lwd=1.6, col=HI, lty=2)
  yb <- max(d$y)*1.12
  segments(median(NE$delta)-NEmad, yb, median(NE$delta)+NEmad, yb, col=HI, lwd=2.2)
  segments(c(median(NE$delta)-NEmad, median(NE$delta)+NEmad), yb-max(d$y)*0.03,
           c(median(NE$delta)-NEmad, median(NE$delta)+NEmad), yb+max(d$y)*0.03,
           col=HI, lwd=2.2)
  text(median(NE$delta), yb+max(d$y)*0.10, sprintf("\u00b1 %.2f", NEmad),
       adj=0.5, cex=1.0, col=HI, font=2)
  axis(1, at=seq(-1.5,1.5,0.5), cex.axis=0.9, lwd=0.6)
  mtext(expression(Delta), side=1, line=2.6, cex=1.05)
  mtext(expression(bold("Nepenthes: true " * Delta * " is 0")), side=3, adj=0,
        line=2.4, cex=1.15)
  mtext(sprintf("%s genes | centre %+.3f (p = %.2f) | sd %.2f, robust sd %.2f",
        format(nrow(NE), big.mark=","), median(NE$delta),
        suppressWarnings(wilcox.test(NE$delta)$p.value), NEsd, NEmad),
        side=3, adj=0, line=0.8, cex=0.82, col=N_COL)

  par(mar=c(4.8,6.4,4.2,1.4))
  SUM <- W[, .(med=median(s)), by=sp][order(med)]
  xl <- c(0, max(c(W$s, PRED))*1.06)
  plot(NA, xlim=xl, ylim=c(0.4, nrow(SUM)+0.6), axes=FALSE, xlab="", ylab="")
  abline(v=PRED, lty=2, col=HI, lwd=1.5)
  for (i in seq_len(nrow(SUM))) {
    v <- W[sp == SUM$sp[i], s]
    segments(quantile(v,.25), i, quantile(v,.75), i, col=adjustcolor(A_COL,alpha.f=0.45), lwd=7, lend=1)
    points(v, jitter(rep(i, length(v)), amount=0.17), pch=19, cex=0.75,
           col=adjustcolor(A_DK, alpha.f=0.55))
    points(median(v), i, pch=18, cex=2.1, col="white")
    points(median(v), i, pch=18, cex=1.6, col=A_DK)
    mtext(SUM$sp[i], side=2, at=i, las=1, line=0.5, cex=0.98, font=3)
    mtext(sprintf("%d units", length(v)), side=4, at=i, las=1, line=-1.4,
          cex=0.75, col=N_COL) }
  axis(1, cex.axis=0.9, lwd=0.6)
  mtext(expression("spread of " * Delta * " inside one block"), side=1, line=2.7, cex=1.0)
  mtext("Measured, not assumed", side=3, adj=0, line=2.4, cex=1.15, font=2)
  mtext(sprintf("one point per species x region x chromosome (%d units, %s genes)",
        nrow(W), format(sum(W$n), big.mark=",")),
        side=3, adj=0, line=0.8, cex=0.82, col=N_COL)
  text(PRED, nrow(SUM)+0.48, " predicted 0.32", adj=0, cex=0.85, col=HI) }
W2 <- 13; H <- 5.2
pdf(file.path(OUT,"PFnoise.pdf"), width=W2, height=H, useDingbats=FALSE); draw(); dev.off()
png(file.path(OUT,"PFnoise.png"), width=W2*160, height=H*160, res=160); draw(); dev.off()
cat("  wrote PFnoise.pdf / .png\n")
