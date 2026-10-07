#!/usr/bin/env Rscript
# ============================================================================
# PF_edges -- are block boundaries real edges, or arbitrary cuts?
#
# THE QUESTION
#   A block is an "atom" of ancestry only if its EDGES are real: Delta should
#   shift abruptly where a boundary is called and stay flat inside a block.
#   If the blocks were arbitrary partitions of a smear, a called boundary would
#   look no different from a random position in a block's interior.
#
# THE TEST
#   At each called A/B boundary take W measured genes either side and record
#   |mean(left) - mean(right)|. Do the same at random INTERIOR positions, at
#   least W genes from any boundary. Aggregate over all boundaries genome-wide.
#     boundary jumps >> interior jumps  -> the edges are real
#     the two overlap                   -> the blocks are arbitrary cuts
#
#   Repeated across W to give a RESOLUTION: the smallest W at which the jump is
#   still separable is how precisely a boundary can be placed.
#
#   Uses raw per-gene Delta only. The boundary positions come from DR02, but
#   the interior comparison controls for that -- both draw on the same blocks.
#
# IN  DR/out/DR02_gene_delta.csv, DR02_propagated_blocks.csv
# OUT presentation_figures/PFedges.{pdf,png}, PFedges_stats.csv
# ============================================================================
suppressMessages(library(data.table))
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
set.seed(23)
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }
A_COL <- "#1D9E75"; A_DK <- "#0F6E56"; B_COL <- "#D85A30"; B_DK <- "#993C1D"
NUL <- "#9C9A92"; N_COL <- "#6E6D69"; HI <- "#BA7517"
WTRY <- c(5, 10, 20, 40); WMAIN <- 20; NINT <- 4000

hr("0. MEASURED GENES, WITH THEIR BLOCK")
G <- fread("DR/out/DR02_gene_delta.csv")[is.finite(delta)]
P <- fread("DR/out/DR02_propagated_blocks.csv")[label %in% c("A","B")]
gk <- if ("gene" %in% names(P)) "gene" else "tip"
M <- merge(G[, .(genome, chr, gene, mid, delta)],
           P[, .(genome, chr, gene=get(gk), region, label)], by=c("genome","chr","gene"))
M[, mb := mid/1e6]; M[, ck := paste(genome, chr)]
setorder(M, ck, mb)
M[, blk := rleid(paste(region,label)), by=ck]
M[, i := seq_len(.N), by=ck]
cat(sprintf("  measured genes with a block: %s | chromosomes %d | blocks %d\n",
    format(nrow(M), big.mark=","), uniqueN(M$ck), nrow(unique(M[, .(ck, blk)]))))
bs <- M[, .N, by=.(ck, blk)]
cat(sprintf("  measured genes per block: median %d | IQR %d-%d\n",
    as.integer(median(bs$N)), as.integer(quantile(bs$N,.25)), as.integer(quantile(bs$N,.75))))

hr("1. JUMP AT BOUNDARIES vs INSIDE BLOCKS")
jump <- function(v, p, w){
  if (p-w < 1 || p+w > length(v)) return(NA_real_)
  abs(mean(v[(p-w+1):p]) - mean(v[(p+1):(p+w)])) }
RES <- rbindlist(lapply(WTRY, function(w){
  bnd <- rbindlist(lapply(split(M, M$ck), function(z){
    z <- z[order(i)]; e <- which(head(z$blk,-1) != tail(z$blk,-1))
    if (!length(e)) return(NULL)
    lab <- vapply(e, function(p) as.integer(z$label[p] != z$label[p+1]), integer(1))
    data.table(ck=z$ck[1], p=e, switch=lab,
               j=vapply(e, function(p) jump(z$delta, p, w), numeric(1))) }))
  bnd <- bnd[is.finite(j)]
  itr <- rbindlist(lapply(split(M, M$ck), function(z){
    z <- z[order(i)]; e <- which(head(z$blk,-1) != tail(z$blk,-1))
    ok <- setdiff((w+1):(nrow(z)-w), unlist(lapply(e, function(p) (p-w):(p+w))))
    if (length(ok) < 2) return(NULL)
    s <- sample(ok, min(length(ok), ceiling(NINT/uniqueN(M$ck))))
    data.table(ck=z$ck[1], j=vapply(s, function(p) jump(z$delta, p, w), numeric(1))) }))
  itr <- itr[is.finite(j)]
  sw <- bnd[switch==1]
  data.table(w=w, n_bnd=nrow(bnd), n_switch=nrow(sw), n_int=nrow(itr),
             bnd_med=median(bnd$j), switch_med=median(sw$j), int_med=median(itr$j),
             ratio=median(sw$j)/median(itr$j),
             p=suppressWarnings(wilcox.test(sw$j, itr$j)$p.value),
             auc={ r <- rank(c(sw$j, itr$j))
                   (sum(r[seq_len(nrow(sw))]) - nrow(sw)*(nrow(sw)+1)/2) /
                     (nrow(sw)*nrow(itr)) }) }))
print(RES[, .(w, n_switch, n_int, switch_med=round(switch_med,3),
              int_med=round(int_med,3), ratio=round(ratio,2),
              auc=round(auc,3), p=signif(p,3))], row.names=FALSE)
cat("\n  switch_med = jump at boundaries where the SUBGENOME changes\n")
cat("  int_med    = jump at random positions inside a block\n")
cat("  auc        = chance a random boundary jumps more than a random interior\n")
cat("               point. 0.5 = indistinguishable, 1.0 = perfectly separated.\n")
fwrite(RES, file.path(OUT,"PFedges_stats.csv"))
best <- RES[auc == max(auc)][1]
cat(sprintf("\n  best separation at W = %d: AUC %.3f, jump %.2fx interior, p = %.3g\n",
    best$w, best$auc, best$ratio, best$p))
if (max(RES$auc) < 0.65)
  cat("  >>> BOUNDARIES ARE NOT DISTINGUISHABLE FROM BLOCK INTERIORS\n")
small <- RES[auc > 0.75][1]
cat(sprintf("  resolution: separable (AUC > 0.75) down to W = %s genes a side\n",
    if (nrow(small) && !is.na(small$w)) as.character(small$w) else "no tested W"))

hr("2. PROFILE ACROSS A BOUNDARY  (subgenome-switching edges only)")
OFF <- -40:40
prof <- rbindlist(lapply(split(M, M$ck), function(z){
  z <- z[order(i)]; e <- which(head(z$blk,-1) != tail(z$blk,-1))
  e <- e[vapply(e, function(p) z$label[p] != z$label[p+1], logical(1))]
  if (!length(e)) return(NULL)
  rbindlist(lapply(e, function(p){
    k <- p + OFF; ok <- k >= 1 & k <= nrow(z)
    sgn <- if (z$label[p] == "A") 1 else -1
    data.table(off=OFF[ok], d=z$delta[k[ok]]*sgn,
               same_region = z$region[p] == z$region[p+1]) })) }))
PR <- prof[, .(m = mean(d), se = sd(d)/sqrt(.N), n = .N), by=off][order(off)]
NB <- prof[, .N] %/% length(OFF)
cat(sprintf("  boundaries profiled: %d | genes contributing: %s\n",
    NB, format(nrow(prof), big.mark=",")))
cat(sprintf("  every profiled edge switches subgenome: %s | of which same region: %.0f%%\n",
    TRUE, 100*mean(prof$same_region)))
cat(sprintf("  left of the edge  %+.3f (should be negative = A)\n", PR[off < 0, weighted.mean(m, n)]))
cat(sprintf("  right of the edge %+.3f (should be positive = B)\n", PR[off > 0, weighted.mean(m, n)]))
cat(sprintf("  mean oriented Delta: left of the edge %+.3f | right %+.3f\n",
    PR[off < 0, weighted.mean(m, n)], PR[off > 0, weighted.mean(m, n)]))

hr("3. FIGURE")
draw <- function(){
  par(mar=c(5.0,5.8,4.4,2.0))
  yl <- range(c(PR$m-1.96*PR$se, PR$m+1.96*PR$se))
  yl <- yl + c(-0.06, 0.10)*diff(yl)
  plot(NA, xlim=range(OFF), ylim=yl, axes=FALSE, xlab="", ylab="")
  abline(h=0, col="grey60", lwd=0.8)
  Lf <- PR[off <= 0]; Rt <- PR[off >= 0]
  polygon(c(Lf$off, rev(Lf$off)), c(Lf$m-1.96*Lf$se, rev(Lf$m+1.96*Lf$se)),
          col=adjustcolor(A_COL, alpha.f=0.26), border=NA)
  polygon(c(Rt$off, rev(Rt$off)), c(Rt$m-1.96*Rt$se, rev(Rt$m+1.96*Rt$se)),
          col=adjustcolor(B_COL, alpha.f=0.26), border=NA)
  lines(Lf$off, Lf$m, col=A_DK, lwd=2.6)
  lines(Rt$off, Rt$m, col=B_DK, lwd=2.6)
  abline(v=0, col="grey35", lwd=1.3, lty=2)
  axis(1, cex.axis=0.95, lwd=0.6); axis(2, las=1, cex.axis=0.95, lwd=0.6)
  mtext("measured genes from the called boundary", side=1, line=3.0, cex=1.0)
  mtext(expression(Delta * ", oriented so the left block is A"), side=2, line=4.0, cex=1.0)
  mtext(expression(bold(Delta ~ "flips at the called edge")), side=3, adj=0, line=2.5, cex=1.25)
  mtext(sprintf("%d boundaries where the subgenome changes, stacked | shading = 95%% interval",
        NB), side=3, adj=0, line=1.0, cex=0.88, col=N_COL)
  text(-38, yl[2]*0.88, "A block", adj=0, cex=1.05, font=2, col=A_DK)
  text(38, yl[1]*0.88, "B block", adj=1, cex=1.05, font=2, col=B_DK)
  text(0, yl[2]*0.97, "called edge", adj=0.5, cex=0.88, col="grey35") }
Wd <- 8.6; H <- 5.4
pdf(file.path(OUT,"PFedges.pdf"), width=Wd, height=H, useDingbats=FALSE); draw(); dev.off()
png(file.path(OUT,"PFedges.png"), width=Wd*160, height=H*160, res=160); draw(); dev.off()
cat("  wrote PFedges.pdf / .png\n")
