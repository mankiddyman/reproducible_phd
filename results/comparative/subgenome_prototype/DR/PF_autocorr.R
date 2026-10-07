#!/usr/bin/env Rscript
# ============================================================================
# PF_autocorr -- does a gene's Delta predict its neighbour's?
#
# THE TEST
#   Correlate Delta at gene i with Delta at gene i+k, along each chromosome,
#   pooled genome-wide. No windows, no thresholds, no block boundaries, no
#   labels -- so nothing here can be circular with DR02's calls.
#     pure noise      -> correlation 0 at every lag
#     real blocks     -> positive at short lag, decaying to 0 at the lag where
#                        genes stop sharing a block
#   The DECAY DISTANCE is an independent estimate of block size.
#
#   Null: Delta shuffled WITHIN each chromosome. Keeps every value and the
#   A:B balance, destroys only the ordering.
#
# EXPECTED SIZE  true means at -0.16 / +0.16 give signal variance ~0.0256;
#   observed Var(Delta) is printed below. rho(1) ~ Var(signal)/Var(total),
#   roughly 0.09-0.20. Small, but with ~13,600 genes the SE is ~0.009.
#
# CAVEAT  only genes with a measured Delta enter, so a lag of 1 means "the next
#   MEASURED gene", not the physically adjacent one. Lags are therefore also
#   reported in Mb so the decay distance can be read physically.
#
# IN  DR/out/DR02_gene_delta.csv, DR02_propagated_blocks.csv (block sizes only)
# OUT presentation_figures/PFautocorr.{pdf,png}, PFautocorr_stats.csv
# ============================================================================
suppressMessages(library(data.table))
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
set.seed(13)
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }
A_COL <- "#1D9E75"; A_DK <- "#0F6E56"; NUL <- "#9C9A92"; N_COL <- "#6E6D69"; HI <- "#BA7517"
NPERM <- 200
LAGS  <- c(1,2,3,5,8,12,20,30,50,80,120,200)

hr("0. INPUT")
G <- fread("DR/out/DR02_gene_delta.csv")[is.finite(delta)]
G[, mb := mid/1e6]; G[, ck := paste(genome, chr)]
setorder(G, ck, mb)
G[, idx := seq_len(.N), by=ck]
n_ck <- G[, .N, by=ck]
cat(sprintf("  measured genes %s | chromosomes %d | per chromosome: median %d (%d-%d)\n",
    format(nrow(G), big.mark=","), nrow(n_ck), as.integer(median(n_ck$N)),
    min(n_ck$N), max(n_ck$N)))
cat(sprintf("  Var(Delta) = %.4f | predicted signal variance 0.16^2 = %.4f\n",
    var(G$delta), 0.16^2))
cat(sprintf("  so rho(1) should be about %.3f if blocks are real\n", 0.16^2/var(G$delta)))

hr("1. AUTOCORRELATION BY LAG")
lagcor <- function(d, k){
  num <- 0; den1 <- 0; den2 <- 0; np <- 0
  for (v in d) { n <- length(v); if (n <= k) next
    a <- v[1:(n-k)]; b <- v[(k+1):n]
    num <- num + sum((a-mean(v))*(b-mean(v)))
    den1 <- den1 + sum((a-mean(v))^2); den2 <- den2 + sum((b-mean(v))^2)
    np <- np + (n-k) }
  if (np < 50) return(c(NA_real_, np))
  c(num/sqrt(den1*den2), np) }
dl <- split(G$delta, G$ck)
AC <- rbindlist(lapply(LAGS, function(k){
  r <- lagcor(dl, k)
  nl <- vapply(seq_len(NPERM), function(b)
    lagcor(lapply(dl, sample), k)[1], numeric(1))
  nl <- nl[is.finite(nl)]
  data.table(lag=k, rho=r[1], pairs=r[2],
             null_med=median(nl), null_lo=quantile(nl,.025), null_hi=quantile(nl,.975),
             z=(r[1]-median(nl))/sd(nl),
             p=(1+sum(abs(nl-median(nl)) >= abs(r[1]-median(nl))))/(length(nl)+1)) }))
print(AC[, .(lag, pairs, rho=round(rho,4), null_hi=round(null_hi,4),
             z=round(z,1), p=signif(p,3))], row.names=FALSE)
fwrite(AC, file.path(OUT,"PFautocorr_stats.csv"))
cat(sprintf("\n  rho at lag 1: %.4f (null 95%% ceiling %.4f, z = %.1f)\n",
    AC$rho[1], AC$null_hi[1], AC$z[1]))
above <- AC[rho > null_hi]
cat(sprintf("  lags above the null envelope: %s\n",
    if (nrow(above)) paste(above$lag, collapse=", ") else "NONE"))
if (!nrow(above)) cat("  >>> no spatial structure detected -- the phasing would not survive this\n")
half <- AC[rho < AC$rho[1]/2][1]
cat(sprintf("  rho falls below half its lag-1 value at lag %s\n",
    if (nrow(half) && !is.na(half$lag)) as.character(half$lag) else "never within the tested range"))

hr("2. DECAY DISTANCE IN PHYSICAL UNITS")
sp <- G[, .(span = diff(range(mb)), n = .N), by=ck][n >= 20]
mb_per <- median(sp$span/sp$n)
cat(sprintf("  median Mb between consecutive MEASURED genes: %.3f\n", mb_per))
for (k in c(1, 5, 20, 50)) cat(sprintf("    lag %3d = about %6.1f Mb\n", k, k*mb_per))
B <- fread("DR/out/DR02_propagated_blocks.csv")[label %in% c("A","B") & is.finite(mb)]
bl <- B[order(mb), { r <- rle(paste(region,label)); en <- cumsum(r$lengths)
                     st <- en - r$lengths + 1
                     .(span = mb[en] - mb[st]) }, by=.(genome, chr)]
cat(sprintf("\n  DR02 block spans: median %.1f Mb | IQR %.1f-%.1f\n",
    median(bl$span), quantile(bl$span,.25), quantile(bl$span,.75)))
cat("  ^ if the autocorrelation decays over a comparable distance, two\n")
cat("    independent routes agree on block size.\n")

hr("3. FIGURE")
draw <- function(){
  par(mar=c(4.8,5.6,4.2,1.8))
  ymax <- max(c(AC$rho, AC$null_hi), na.rm=TRUE)*1.20
  ymin <- min(c(0, AC$null_lo, AC$rho), na.rm=TRUE)*1.5
  plot(NA, xlim=c(min(LAGS), max(LAGS)), ylim=c(ymin, ymax), log="x",
       axes=FALSE, xlab="", ylab="")
  polygon(c(AC$lag, rev(AC$lag)), c(AC$null_hi, rev(AC$null_lo)),
          col=adjustcolor(NUL, alpha.f=0.40), border=NA)
  abline(h=0, col="grey55", lwd=0.8)
  lines(AC$lag, AC$rho, col=A_DK, lwd=2.6)
  points(AC$lag, AC$rho, pch=19, cex=1.4,
         col=ifelse(AC$rho > AC$null_hi, A_COL, NUL))
  axis(1, at=LAGS, labels=LAGS, cex.axis=0.85, lwd=0.6)
  axis(2, las=1, cex.axis=0.9, lwd=0.6)
  mtext("distance between genes (measured genes apart)", side=1, line=2.9, cex=0.98)
  mtext("correlation between their \u0394", side=2, line=4.0, cex=0.98)
  mtext("A gene's vote predicts its neighbour's", side=3, adj=0, line=2.4,
        cex=1.25, font=2)
  mtext(sprintf("%s genes, every chromosome | grey band = \u0394 shuffled within each chromosome, %d times",
        format(nrow(G), big.mark=","), NPERM),
        side=3, adj=0, line=0.9, cex=0.88, col=N_COL)
  text(LAGS[1]*1.1, AC$rho[1], sprintf("  %.3f", AC$rho[1]),
       adj=0, cex=1.0, font=2, col=A_DK)
  zl <- AC[rho > null_hi]
  if (nrow(zl)) text(max(LAGS)*0.95, ymax*0.88,
    sprintf("above chance out to %d genes apart\n(about %.0f Mb)",
            max(zl$lag), max(zl$lag)*mb_per), adj=1, cex=0.95, col=HI) }
W <- 8.8; H <- 5.4
pdf(file.path(OUT,"PFautocorr.pdf"), width=W, height=H, useDingbats=FALSE); draw(); dev.off()
png(file.path(OUT,"PFautocorr.png"), width=W*160, height=H*160, res=160); draw(); dev.off()
cat("  wrote PFautocorr.pdf / .png\n")
