#!/usr/bin/env Rscript
# ============================================================================
# PF_phasing -- how per-gene Delta votes become phased blocks
#
# LEFT, one chromosome, three tracks:
#   1 raw Delta for every gene, coloured by its own vote. Roughly a quarter sit
#     on the wrong side of zero -- that is the predicted noise, not an error.
#   2 the same genes averaged per block: medians with 95% intervals, far apart.
#   3 the resulting A/B blocks, with ancestral-region boundaries marked in
#     GENESPACE's own region colours (taken from the saved riparian object).
#
# RIGHT how often a window of n genes reproduces its block's call:
#   1 gene 75%, 10 genes 92%, 25 genes 96%, 100 genes 98%.
#   NOTE this is internal consistency -- the window is compared against a call
#   derived from the same genes -- so it measures stability, not accuracy.
#   Windows are confined WITHIN a block, so they never straddle a boundary.
#
# Blocks are run-length encodings of (ancestral region, A/B label) along the
# chromosome, gene by gene. DR02_segments.csv is NOT used: its mb_lo/mb_hi are
# per-region min/max, so a region occurring in several blocks gets one box
# spanning everything between them.
#
# IN  DR/out/DR02_gene_delta.csv, DR02_propagated_blocks.csv,
#     ../genespace/riparian/Nepenthes_gracilis_geneOrder_rSourceData.rda
# OUT presentation_figures/PFphasing.{pdf,png}, PFphasing_accuracy.csv
# ============================================================================
suppressMessages(library(data.table))
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
set.seed(3)
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
GS <- file.path(dirname(getwd()), "genespace")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }

A_COL <- "#1D9E75"; B_COL <- "#D85A30"; A_DK <- "#0F6E56"; B_DK <- "#993C1D"
N_COL <- "#6E6D69"
EX <- c("Drosera_paradoxa", "chr4_hap1")
WSZ <- c(1, 5, 10, 25, 50, 100)

hr("0. BLOCKS, GENE BY GENE")
P <- fread("DR/out/DR02_propagated_blocks.csv")[label %in% c("A","B") & is.finite(mb)]
gk <- if ("gene" %in% names(P)) "gene" else "tip"
blocks <- function(d){ d <- d[order(mb)]
  r <- rle(paste(d$region, d$label)); en <- cumsum(r$lengths); st <- en - r$lengths + 1
  data.table(region=sub(" .*$","",r$values), label=sub("^.* ","",r$values),
             n=r$lengths, lo=d$mb[st], hi=d$mb[en]) }
TR <- P[, blocks(.SD), by=.(genome, chr), .SDcols=c("mb","region","label")]
cat(sprintf("  genes %s | blocks %d | singletons %d | all >= 10 genes: %s\n",
    format(nrow(P), big.mark=","), nrow(TR), sum(TR$n==1), all(TR$n>=10)))

G <- fread("DR/out/DR02_gene_delta.csv")[is.finite(delta)]; G[, mb := mid/1e6]
M <- merge(G[, .(genome, chr, gene, mb, delta)],
           P[, .(genome, chr, gene=get(gk), region, label)], by=c("genome","chr","gene"))
M[, vote := fifelse(delta<0, "A", "B")]
cat(sprintf("  genes with a measured Delta: %s | agreeing with their block: %.1f%%\n",
    format(nrow(M), big.mark=","), 100*mean(M$vote==M$label)))
cat(sprintf("  predicted disagreement from the noise model: %.1f%% | observed %.1f%%\n",
    100*pnorm(0.5, lower.tail=FALSE), 100*mean(M$vote!=M$label)))

hr("1. ACCURACY BY WINDOW SIZE  (within blocks only)")
M[, blk := paste(genome, chr, region, label)]
AC <- rbindlist(lapply(WSZ, function(w){
  r <- M[, { v <- .SD[order(mb)]
    if (.N < w) NULL else { k <- floor(.N/w)
      g <- rep(seq_len(k), each=w)[seq_len(k*w)]
      d <- data.table(dl=v$delta[seq_len(k*w)], lb=v$label[seq_len(k*w)], g=g)
      d[, .(ok = (median(dl)<0) == (lb[1]=="A")), by=g] } }, by=blk]
  data.table(w=w, windows=nrow(r), acc=mean(r$ok)) }))
print(AC, row.names=FALSE)
if (is.unsorted(AC$acc)) cat("  !! accuracy not monotone -- check the windowing\n")
fwrite(AC, file.path(OUT,"PFphasing_accuracy.csv"))

hr("2. GENESPACE REGION COLOURS")
RDA <- file.path(GS,"riparian","Nepenthes_gracilis_geneOrder_rSourceData.rda")
RCOL <- NULL
if (file.exists(RDA)) {
  e <- new.env(); load(RDA, envir=e)
  BR <- as.data.frame(e$srcd$ggplotObj$layers[[1]]$data)
  BR$region <- sub("^.*_[0-9]+ ", "", BR$blkID)
  rc <- unique(BR[, c("region","color")])
  rc <- rc[!duplicated(rc$region), ]
  RCOL <- setNames(rc$color, rc$region)
  cat("  regions with a GENESPACE colour:", length(RCOL), "\n"); print(RCOL)
} else cat("  riparian object not found -- falling back to grey region ticks\n")

hr("3. FIGURE")
E  <- P[genome==EX[1] & chr==EX[2]][order(mb)]
ET <- TR[genome==EX[1] & chr==EX[2]][order(lo)]
GE <- M[genome==EX[1] & chr==EX[2]]
if (!nrow(ET)) die("no blocks for %s %s", EX[1], EX[2])
cat(sprintf("  %s %s: %s genes, %d blocks, %d regions\n",
    sub("Drosera_","D. ",EX[1]), EX[2], format(nrow(E), big.mark=","),
    nrow(ET), uniqueN(ET$region)))
ST <- rbindlist(lapply(seq_len(nrow(ET)), function(i){
  v <- GE[mb>=ET$lo[i] & mb<=ET$hi[i], delta]
  if (length(v) < 5) return(NULL)
  data.table(i=i, m=median(v), se=mad(v)/sqrt(length(v)), n=length(v),
             lo=ET$lo[i], hi=ET$hi[i], label=ET$label[i], region=ET$region[i]) }))
cat(sprintf("  block medians: A %.2f..%.2f | B %.2f..%.2f | gap %.2f\n",
    min(ST$m[ST$label=="A"]), max(ST$m[ST$label=="A"]),
    min(ST$m[ST$label=="B"]), max(ST$m[ST$label=="B"]),
    min(ST$m[ST$label=="B"]) - max(ST$m[ST$label=="A"])))

left <- function(){
  xr <- range(E$mb); par(mar=c(4.2,1.4,3.4,1.6))
  plot(NA, xlim=xr, ylim=c(-1.30, 3.30), axes=FALSE, xlab="", ylab="")

  segments(xr[1], 2.30, xr[2], 2.30, col="grey70", lwd=0.7, lty=2)
  points(GE$mb, 2.30 + GE$delta*0.72, pch=19, cex=0.30,
         col=adjustcolor(ifelse(GE$delta<0, A_COL, B_COL), alpha.f=0.30))
  text(xr[1], 3.22, "1  every gene votes", adj=0, cex=1.05, font=2)
  text(xr[1], 3.02, sprintf("%.0f%% land on the wrong side \u2014 that is the noise",
       100*mean(GE$vote!=GE$label)), adj=0, cex=0.86, col=N_COL)
  text(xr[2], 2.30, " 0", adj=0, cex=0.8, col=N_COL, xpd=NA)

  rng <- max(abs(c(ST$m-1.96*ST$se, ST$m+1.96*ST$se)))
  yy <- function(v) 1.00 + v/rng*0.62
  segments(xr[1], yy(0), xr[2], yy(0), col="grey70", lwd=0.7, lty=2)
  for (k in seq_len(nrow(ST))) {
    cl <- if (ST$label[k]=="A") A_DK else B_DK
    rect(ST$lo[k], yy(ST$m[k]-1.96*ST$se[k]), ST$hi[k], yy(ST$m[k]+1.96*ST$se[k]),
         col=adjustcolor(cl, alpha.f=0.28), border=NA)
    segments(ST$lo[k], yy(ST$m[k]), ST$hi[k], yy(ST$m[k]), col=cl, lwd=4.5, lend=1) }
  text(xr[1], 1.74, "2  averaged over a block, the two are far apart",
       adj=0, cex=1.05, font=2)
  text(xr[2], yy(0), " 0", adj=0, cex=0.8, col=N_COL, xpd=NA)

  for (i in seq_len(nrow(ET)))
    rect(ET$lo[i], -0.60, ET$hi[i], -0.16,
         col=if (ET$label[i]=="A") A_COL else B_COL, border="white", lwd=0.6)
  for (i in seq_len(nrow(ET)))
    if (ET$hi[i]-ET$lo[i] > diff(xr)*0.035)
      text(mean(c(ET$lo[i],ET$hi[i])), -0.38, ET$label[i],
           adj=0.5, cex=1.0, font=2, col="white")
  for (i in seq_len(nrow(ET))) {
    cl <- if (!is.null(RCOL) && ET$region[i] %in% names(RCOL)) RCOL[[ET$region[i]]] else "grey70"
    rect(ET$lo[i], -0.78, ET$hi[i], -0.64, col=cl, border="white", lwd=0.6) }
  text(xr[1], 0.00, "3  each block is called", adj=0, cex=1.05, font=2)
  text(xr[2], -0.71, " ancestral", adj=0, cex=0.72, col=N_COL, xpd=NA)

  axis(1, cex.axis=0.9, lwd=0.6, pos=-1.00)
  mtext("position (Mb)", side=1, line=2.6, cex=0.9)
  mtext(sprintf("%s  %s", sub("Drosera_","D. ",EX[1]), EX[2]),
        side=3, adj=0, line=1.6, cex=1.2, font=2)
  mtext(sprintf("%s genes | %d blocks | a block ends where the subgenome or the ancestral chromosome changes",
        format(nrow(E), big.mark=","), nrow(ET)),
        side=3, adj=0, line=0.4, cex=0.84, col=N_COL) }

right <- function(){
  par(mar=c(4.6,5.2,3.4,1.6))
  plot(NA, xlim=c(0.8, max(WSZ)*1.25), ylim=c(70, 100), log="x", axes=FALSE,
       xlab="", ylab="")
  abline(h=seq(70,100,5), col="grey92", lwd=0.6)
  lines(AC$w, 100*AC$acc, col=A_DK, lwd=2.4)
  points(AC$w, 100*AC$acc, pch=19, cex=1.5, col=A_COL)
  for (i in seq_len(nrow(AC)))
    text(AC$w[i], 100*AC$acc[i] + 2.0, sprintf("%.0f%%", 100*AC$acc[i]),
         adj=0.5, cex=0.92, col=A_DK)
  axis(1, at=WSZ, labels=WSZ, cex.axis=0.9, lwd=0.6)
  axis(2, at=seq(70,100,10), labels=paste0(seq(70,100,10),"%"), las=1,
       cex.axis=0.9, lwd=0.6)
  mtext("genes averaged", side=1, line=2.7, cex=0.95)
  mtext("agreement with the block call", side=2, line=3.6, cex=0.95)
  mtext("One gene is a guess. A block is not.", side=3, adj=0, line=1.6,
        cex=1.2, font=2)
  mtext(sprintf("%s genes across %d blocks | windows never cross a block boundary",
        format(nrow(M), big.mark=","), uniqueN(M$blk)),
        side=3, adj=0, line=0.4, cex=0.84, col=N_COL) }

draw <- function(){ layout(matrix(1:2, nrow=1), widths=c(1.5,1)); left(); right() }
W <- 14; H <- 5.6
pdf(file.path(OUT,"PFphasing.pdf"), width=W, height=H, useDingbats=FALSE); draw(); dev.off()
png(file.path(OUT,"PFphasing.png"), width=W*160, height=H*160, res=160); draw(); dev.off()
cat("  wrote PFphasing.pdf / .png\n")
