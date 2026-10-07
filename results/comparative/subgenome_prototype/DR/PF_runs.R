#!/usr/bin/env Rscript
# ============================================================================
# PF_runs -- does Delta hold one sign over long stretches, more than chance?
#
# THE TEST, and why it is not circular:
#   One caller, applied identically to real and shuffled data:
#     tile the chromosome into non-overlapping windows of W genes
#     call each window by the SIGN of its median raw Delta
#     run-length encode the window calls
#   It never reads a DR02 label, never uses a DR02 boundary, and does not score
#   genes against a call derived from those same genes. The ONLY difference
#   between real and null is whether Delta sits where it actually sits:
#   the null shuffles Delta among the genes of that chromosome, so every value
#   and the overall A:B balance are preserved and only the ordering is lost.
#
#   Long runs in the real data and short runs in the shuffle = Delta is
#   spatially organised into blocks. Equal runs = the blocks are an artefact.
#
# IN  DR/out/DR02_gene_delta.csv
# OUT presentation_figures/PFruns.{pdf,png}, PFruns_stats.csv
# ============================================================================
suppressMessages(library(data.table))
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
set.seed(11)
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }
A_COL <- "#1D9E75"; B_COL <- "#D85A30"; A_DK <- "#0F6E56"; NUL <- "#9C9A92"; N_COL <- "#6E6D69"
NPERM <- 999
WTRY  <- c(5, 10, 20, 40)
WUSE  <- 10
MINW  <- 8      # a chromosome needs at least this many windows to be testable

hr("0. INPUT")
G <- fread("DR/out/DR02_gene_delta.csv")[is.finite(delta)]
G[, mb := mid/1e6]; G[, ck := paste(genome, chr)]
setorder(G, ck, mb)
cat(sprintf("  genes with a raw Delta: %s across %d chromosomes\n",
    format(nrow(G), big.mark=","), uniqueN(G$ck)))
pc <- G[, .N, by=ck]
cat(sprintf("  genes per chromosome: median %d | range %d-%d\n",
    as.integer(median(pc$N)), min(pc$N), max(pc$N)))

# one caller, used for real and shuffled alike
callruns <- function(d, w){
  k <- floor(length(d)/w)
  if (k < 2) return(NULL)
  m <- vapply(seq_len(k), function(i) median(d[((i-1)*w+1):(i*w)]), numeric(1))
  s <- ifelse(m < 0, "A", "B")
  list(calls=s, nwin=k, longest=max(rle(s)$lengths), nruns=length(rle(s)$lengths)) }

hr("1. WINDOW SIZE -- which has enough windows to test?")
for (w in WTRY) {
  ok <- pc[floor(N/w) >= MINW]
  cat(sprintf("  W = %2d genes | chromosomes with >= %d windows: %2d of %d | median windows %.0f\n",
      w, MINW, nrow(ok), nrow(pc),
      if (nrow(ok)) median(floor(ok$N/w)) else 0)) }
cat(sprintf("\n  using W = %d\n", WUSE))

hr("2. RUN-LENGTH, REAL vs SHUFFLED")
cks <- pc[floor(N/WUSE) >= MINW, ck]
if (!length(cks)) die("no chromosome has enough windows at W = %d", WUSE)
R <- rbindlist(lapply(cks, function(k){
  d <- G[ck==k, delta]
  o <- callruns(d, WUSE); if (is.null(o)) return(NULL)
  nl <- vapply(seq_len(NPERM), function(b){
    z <- callruns(sample(d), WUSE); if (is.null(z)) NA_real_ else z$longest }, numeric(1))
  nl <- nl[is.finite(nl)]
  data.table(chrom=k, genes=length(d), windows=o$nwin,
             obs=o$longest, nruns=o$nruns,
             null_med=median(nl), null_p95=as.numeric(quantile(nl,.95)),
             p=(1+sum(nl >= o$longest))/(length(nl)+1)) }))
R[, ratio := obs/pmax(null_med,1)]
setorder(R, -ratio)
cat(sprintf("  chromosomes tested: %d\n", nrow(R)))
print(head(R, 10), row.names=FALSE)
cat(sprintf("\n  longest run of same-sign windows: real median %.0f | shuffled median %.0f\n",
    median(R$obs), median(R$null_med)))
cat(sprintf("  ratio real:shuffled -- median %.2fx | range %.2fx-%.2fx\n",
    median(R$ratio), min(R$ratio), max(R$ratio)))
cat(sprintf("  real exceeds the shuffled 95th percentile on %d of %d chromosomes\n",
    sum(R$obs > R$null_p95), nrow(R)))
cat(sprintf("  p <= 0.05 on %d | p <= 0.01 on %d\n", sum(R$p<=0.05), sum(R$p<=0.01)))
cat(sprintf("  combined (Fisher) across chromosomes: X2 = %.1f, df = %d, p = %.3g\n",
    -2*sum(log(R$p)), 2*nrow(R),
    pchisq(-2*sum(log(R$p)), 2*nrow(R), lower.tail=FALSE)))
fwrite(R, file.path(OUT,"PFruns_stats.csv"))

hr("3. EXEMPLAR STRIP")
EX <- R$chrom[1]
dex <- G[ck==EX, delta]
oex <- callruns(dex, WUSE)
sex <- lapply(1:3, function(i) callruns(sample(dex), WUSE)$calls)
cat(sprintf("  %s: %d genes, %d windows | real longest run %d | shuffles %s\n",
    EX, length(dex), oex$nwin, oex$longest,
    paste(vapply(sex, function(z) max(rle(z)$lengths), numeric(1)), collapse=", ")))

hr("4. FIGURE")
strip <- function(){
  par(mar=c(3.4,7.0,4.0,1.4))
  n <- oex$nwin
  plot(NA, xlim=c(0.5, n+0.5), ylim=c(-3.9, 0.9), axes=FALSE, xlab="", ylab="")
  dr <- function(v, y, col){
    for (i in seq_along(v))
      rect(i-0.48, y-0.34, i+0.48, y+0.34,
           col=if (v[i]=="A") A_COL else B_COL, border="white", lwd=0.4) }
  dr(oex$calls, 0)
  for (j in 1:3) dr(sex[[j]], -j)
  mtext("real", side=2, at=0, las=1, line=0.6, cex=1.05, font=2)
  for (j in 1:3) mtext(sprintf("shuffle %d", j), side=2, at=-j, las=1, line=0.6,
                       cex=0.95, col=N_COL)
  axis(1, at=pretty(c(1,n)), cex.axis=0.88, lwd=0.6)
  mtext(sprintf("window along the chromosome (%d genes each)", WUSE),
        side=1, line=2.2, cex=0.92)
  mtext("Same genes, same caller, only the order changed",
        side=3, adj=0, line=2.2, cex=1.2, font=2)
  mtext(sprintf("%s | longest same-sign run: %d real, %s shuffled",
        sub("Drosera_","D. ", EX), oex$longest,
        paste(vapply(sex, function(z) max(rle(z)$lengths), numeric(1)), collapse="/")),
        side=3, adj=0, line=0.7, cex=0.88, col=N_COL) }

scat <- function(){
  par(mar=c(4.8,5.4,4.0,1.6))
  lim <- c(1, max(c(R$obs, R$null_med))*1.35)
  plot(NA, xlim=lim, ylim=lim, log="xy", axes=FALSE, xlab="", ylab="")
  abline(0, 1, col="grey60", lwd=1, lty=2)
  points(R$null_med, R$obs, pch=19, cex=1.25,
         col=adjustcolor(ifelse(R$p<=0.05, A_COL, NUL), alpha.f=0.75))
  at <- c(1,2,3,5,10,20,40,80); at <- at[at>=lim[1] & at<=lim[2]]
  axis(1, at=at, labels=at, cex.axis=0.9, lwd=0.6)
  axis(2, at=at, labels=at, las=1, cex.axis=0.9, lwd=0.6)
  mtext("longest run when Delta is shuffled", side=1, line=2.9, cex=0.95)
  mtext("longest run, real", side=2, line=3.6, cex=0.95)
  mtext("Delta is organised into blocks", side=3, adj=0, line=2.2, cex=1.2, font=2)
  mtext(sprintf("%d chromosomes | %d shuffles each | coloured where p <= 0.05",
        nrow(R), NPERM), side=3, adj=0, line=0.7, cex=0.88, col=N_COL)
  text(lim[2]*0.95, lim[1]*1.3, "equal", adj=1, cex=0.88, col="grey45")
  text(lim[1]*1.25, lim[2]*0.80, sprintf("%.1fx longer\n(median)", median(R$ratio)),
       adj=0, cex=1.05, col=A_DK) }

draw <- function(){ layout(matrix(1:2, nrow=1), widths=c(1.3,1)); strip(); scat() }
W <- 13.5; H <- 5.2
pdf(file.path(OUT,"PFruns.pdf"), width=W, height=H, useDingbats=FALSE); draw(); dev.off()
png(file.path(OUT,"PFruns.png"), width=W*160, height=H*160, res=160); draw(); dev.off()
cat("  wrote PFruns.pdf / .png\n")
