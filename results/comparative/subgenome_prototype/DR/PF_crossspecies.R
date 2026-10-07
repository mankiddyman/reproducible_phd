#!/usr/bin/env Rscript
# ============================================================================
# PF_crossspecies -- do A copies resemble A copies, and B resemble B, across
#                    species? Scored SEPARATELY for the two subgenomes.
#
# Delta was computed only from Drosera-to-Dionaea distances, so Drosera-to-
# Drosera distances are free evidence.
#
# THE TEST, per subgenome
#   At a locus where two Drosera species each carry an A copy and a B copy:
#     A test:  is d(s1_A, s2_A) < d(s1_A, s2_B) ?   does A prefer A?
#     B test:  is d(s1_B, s2_B) < d(s1_B, s2_A) ?   does B prefer B?
#   Each is one comparison of two distances sharing a tip, so under random
#   labels each wins half the time. Chance is 50% and AAB does not change that
#   -- the 2:1 composition decides which loci are usable, not the null.
#
#   Scoring them apart is the point: B is the single-copy subgenome and may be
#   phased less well than A. The earlier combined version could not show that.
#
# CONTROL  the same tests with A/B swapped at random within each species, which
#   must land on 50%.
#
# CAVEAT  a locus needs BOTH an A and a B copy in BOTH species, so this runs on
#   a subset of loci. That restricts coverage, not the null.
#
# IN  DR/out/pairwise_ks.csv, DR02_gene_delta.csv, DR02_propagated_blocks.csv
# OUT presentation_figures/PFcross.{pdf,png}, PFcross_stats.csv
# ============================================================================
suppressMessages(library(data.table))
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
set.seed(17)
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }
A_COL <- "#1D9E75"; A_DK <- "#0F6E56"; B_COL <- "#D85A30"; B_DK <- "#993C1D"
NUL <- "#9C9A92"; N_COL <- "#6E6D69"; HI <- "#BA7517"
DSMAX <- 5

hr("0. LABEL EVERY COPY")
G <- fread("DR/out/DR02_gene_delta.csv")[is.finite(delta)]
P <- fread("DR/out/DR02_propagated_blocks.csv")[label %in% c("A","B")]
gk <- if ("gene" %in% names(P)) "gene" else "tip"
L <- unique(merge(G[, .(locus, genome, gene, tip)],
                  P[, .(genome, gene=get(gk), label)],
                  by=c("genome","gene"))[, .(locus, genome, tip, label)])
cat(sprintf("  labelled copies: %s | loci %s\n",
    format(nrow(L), big.mark=","), format(uniqueN(L$locus), big.mark=",")))
ONE <- L[, .(tip = tip[1]), by=.(locus, genome, label)]
W <- dcast(ONE, locus + genome ~ label, value.var="tip")
W <- W[!is.na(A) & !is.na(B)]
setnames(W, c("A","B"), c("tipA","tipB"))
cat(sprintf("  (locus, species) with both an A and a B copy: %s\n",
    format(nrow(W), big.mark=",")))

hr("1. DISTANCES")
K <- fread("DR/out/pairwise_ks.csv")[is.finite(dS) & dS >= 0 & dS < DSMAX & codons >= 100]
K <- K[!(sp1 %in% c("Dionaea_muscipula","Nepenthes_gracilis")) &
       !(sp2 %in% c("Dionaea_muscipula","Nepenthes_gracilis"))]
cat(sprintf("  Drosera-Drosera dS rows: %s\n", format(nrow(K), big.mark=",")))
DD <- rbind(K[, .(locus=anchor, a=seq1, b=seq2, d=dS)],
            K[, .(locus=anchor, a=seq2, b=seq1, d=dS)])
dm <- setNames(DD$d, paste(DD$locus, DD$a, DD$b))
gd <- function(l,x,y) unname(dm[paste(l,x,y)])

hr("2. THE TWO TESTS")
run <- function(WD, tag){
  out <- rbindlist(lapply(split(WD, WD$locus), function(z){
    if (nrow(z) < 2) return(NULL)
    cb <- combn(nrow(z), 2)
    rbindlist(lapply(seq_len(ncol(cb)), function(j){
      i1 <- cb[1,j]; i2 <- cb[2,j]; lo <- z$locus[1]
      aa <- gd(lo, z$tipA[i1], z$tipA[i2]); ab <- gd(lo, z$tipA[i1], z$tipB[i2])
      bb <- gd(lo, z$tipB[i1], z$tipB[i2]); ba <- gd(lo, z$tipB[i1], z$tipA[i2])
      rbindlist(list(
        if (all(is.finite(c(aa,ab)))) data.table(locus=lo, sp1=z$genome[i1],
          sp2=z$genome[i2], sub="A", win=aa < ab) else NULL,
        if (all(is.finite(c(bb,ba)))) data.table(locus=lo, sp1=z$genome[i1],
          sp2=z$genome[i2], sub="B", win=bb < ba) else NULL)) })) }))
  if (!nrow(out)) return(NULL)
  out[, pair := paste(pmin(sp1,sp2), pmax(sp1,sp2), sep=" / ")]
  out[, set := tag]; out }

RE <- run(W, "real")
if (is.null(RE)) die("no usable quartets")
WS <- copy(W)
WS[, c("tipA","tipB") := { sw <- runif(.N) < 0.5
  list(fifelse(sw, tipB, tipA), fifelse(sw, tipA, tipB)) }]
SH <- run(WS, "shuffled")

for (sb in c("A","B")) {
  z <- RE[sub==sb]; bt <- binom.test(sum(z$win), nrow(z), 0.5)
  cat(sprintf("  %s tests: %s | same-subgenome closer in %.4f  CI [%.3f, %.3f]  p = %.3g\n",
      sb, format(nrow(z), big.mark=","), mean(z$win),
      bt$conf.int[1], bt$conf.int[2], bt$p.value)) }
if (!is.null(SH)) for (sb in c("A","B")) {
  z <- SH[sub==sb]
  cat(sprintf("  %s CONTROL (labels swapped): %.4f  (must be ~0.500)\n", sb, mean(z$win))) }
dif <- binom.test(sum(RE[sub=="A", win]), nrow(RE[sub=="A"]),
                  mean(RE[sub=="B", win]))$p.value
cat(sprintf("\n  A vs B difference: %.4f vs %.4f | p = %.3g %s\n",
    mean(RE[sub=="A", win]), mean(RE[sub=="B", win]), dif,
    if (dif < 0.05) "-- the two subgenomes are NOT equally well phased" else "-- indistinguishable"))

BP <- RE[, .(tests=.N, frac=mean(win)), by=.(pair, sub)]
BP[, se := sqrt(frac*(1-frac)/tests)]
cat("\n  by species pair:\n")
print(dcast(BP, pair ~ sub, value.var=c("frac","tests"))[order(-frac_A)], row.names=FALSE)
fwrite(rbind(RE, SH, fill=TRUE), file.path(OUT,"PFcross_stats.csv"))

hr("3. FIGURE")
POOL <- RE[, .(frac=mean(win), tests=.N), by=sub]
POOL[, se := sqrt(frac*(1-frac)/tests)]
BP[, grp := fifelse(grepl("regia", pair), "with D. regia", "core Drosera only")]
ordp <- BP[sub=="A"][order(grp, frac), pair]
GRP  <- BP[sub=="A"][order(grp, frac), grp]
cat("\n  grouped:\n")
print(BP[sub=="A", .(pairs=.N, frac=round(mean(frac),3)), by=grp], row.names=FALSE)
draw <- function(){
  par(mar=c(6.4,13.4,5.0,3.4), xpd=NA)
  n <- length(ordp)
  gap <- cumsum(c(0, as.integer(head(GRP,-1) != tail(GRP,-1)))) * 0.9
  ypos <- seq_along(ordp) + gap
  xl <- c(0.45, max(BP$frac + 2*BP$se, 0.96))
  plot(NA, xlim=xl, ylim=c(-0.3, max(ypos)+2.2), axes=FALSE, xlab="", ylab="")
  abline(v=0.5, lty=2, col="grey55", lwd=1.2)
  for (g in unique(GRP)) {
    yy <- ypos[GRP == g]
    rect(xl[1], min(yy)-0.45, xl[2], max(yy)+0.45,
         col=adjustcolor(if (grepl("regia", g)) HI else N_COL, alpha.f=0.05), border=NA)
    text(xl[1]-0.005, max(yy)+0.42, g, adj=c(1,0), cex=0.92, font=3,
         col=if (grepl("regia", g)) HI else N_COL) }
  for (i in seq_along(ordp)) {
    for (sb in c("A","B")) {
      z <- BP[pair==ordp[i] & sub==sb]; if (!nrow(z)) next
      off <- if (sb=="A") 0.17 else -0.17
      cl  <- if (sb=="A") A_COL else B_COL
      segments(z$frac-1.96*z$se, ypos[i]+off, z$frac+1.96*z$se, ypos[i]+off,
               col=adjustcolor(cl, alpha.f=0.45), lwd=2.6, lend=1)
      points(z$frac, ypos[i]+off, pch=19, cex=1.55, col=cl) }
    mtext(sub("Drosera_","D. ", sub(" / Drosera_", " / D. ", ordp[i])),
          side=2, at=ypos[i], las=1, line=0.6, cex=0.95)
    mtext(format(BP[pair==ordp[i] & sub=="A", tests], big.mark=","), side=4,
          at=ypos[i], las=1, line=0.3, cex=0.78, col=N_COL) }
  ya <- max(ypos) + 1.5
  for (sb in c("A","B")) {
    z <- POOL[sub==sb]; off <- if (sb=="A") 0.22 else -0.22
    cl <- if (sb=="A") A_DK else B_DK
    segments(z$frac-1.96*z$se, ya+off, z$frac+1.96*z$se, ya+off, col=cl, lwd=3.4, lend=1)
    points(z$frac, ya+off, pch=18, cex=2.4, col=cl)
    text(z$frac+0.012, ya+off, sprintf("%.0f%%", 100*z$frac),
         adj=0, cex=1.0, font=2, col=cl) }
  mtext("all pairs", side=2, at=ya, las=1, line=0.6, cex=1.05, font=2)
  axis(1, at=seq(0.5,0.95,0.1), labels=sprintf("%.0f%%", 100*seq(0.5,0.95,0.1)),
       cex.axis=0.95, lwd=0.6, pos=-0.25)
  mtext("same-subgenome copies are the closer pair", side=1, line=3.4, cex=1.0)
  text(0.5, -0.9, "chance", adj=0.5, cex=0.88, col="grey45")
  legend(xl[1], -1.7, bty="n", cex=1.0, horiz=TRUE, xpd=NA,
         legend=c("A subgenome","B subgenome"), pch=19, col=c(A_COL,B_COL))
  mtext("A resembles A, B resembles B", side=3, adj=0, line=3.1, cex=1.25, font=2)
  mtext(sprintf("%s tests | Drosera-to-Drosera distances, never used to make the labels",
        format(nrow(RE), big.mark=",")), side=3, adj=0, line=1.7, cex=0.88, col=N_COL)
  if (!is.null(SH))
    mtext(sprintf("labels swapped: A %.0f%%, B %.0f%%",
          100*mean(SH[sub=="A", win]), 100*mean(SH[sub=="B", win])),
          side=3, adj=0, line=0.4, cex=0.86, col=HI) }
Wd <- 10.8; H <- 6.4
pdf(file.path(OUT,"PFcross.pdf"), width=Wd, height=H, useDingbats=FALSE); draw(); dev.off()
png(file.path(OUT,"PFcross.png"), width=Wd*160, height=H*160, res=160); draw(); dev.off()
cat("  wrote PFcross.pdf / .png\n")
