#!/usr/bin/env Rscript
# ============================================================================
# PF_dionaea -- Dionaea allopolyploidy, one slide: a schematic and a plot
#
#   PFD_1_votes  how each ancestral gene votes (schematic, drawn in R)
#   PFD_2_bands  16 chromosomes fall into two non-overlapping retention bands,
#                each pair beside a grey null dumbbell from the same counts
#   PFD_ALL      the two side by side
#
# retention(member) = (genes it kept) / (genes considered for that pair)
#                   = (k + n_1to2) / (total + n_1to2)
# Green is defined as the member that kept more, so green-above-orange is by
# construction. The GAP is not: under random loss the two members differ only
# by sampling noise, which the grey dumbbells show at the same n.
#
# sg1/sg2 are arbitrary assembly labels -- this shows every pair is asymmetric,
# NOT that one subgenome is dominant genome-wide. That link is an assumption.
#
# IN   fractionation_by_chrpair.csv
# OUT  presentation_figures/PFD_{1_votes,2_bands,ALL}.{pdf,png}
#      presentation_figures/PFD_retention.csv
# ============================================================================
suppressMessages(library(data.table))
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
set.seed(11)
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
hr  <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }

WIN <- "#1D9E75"; LOSE <- "#D85A30"; NUL <- "#9C9A92"; GRY <- "#777770"
emit <- function(nm, w, h, fn) {
  pdf(file.path(OUT, paste0(nm,".pdf")), width=w, height=h, useDingbats=FALSE); fn(); dev.off()
  png(file.path(OUT, paste0(nm,".png")), width=w*160, height=h*160, res=160); fn(); dev.off()
  cat("  wrote ", nm, "\n", sep="") }

hr("0. INPUT + CHECKS")
P <- as.data.frame(fread("fractionation_by_chrpair.csv"))
need <- c("exp_pair","chrA","chrB","kA","kB","total","frac_A","p_adj","n_1to2")
if (length(setdiff(need, names(P)))) die("missing: %s", paste(setdiff(need,names(P)), collapse=", "))
cat("  kA + kB == total   : ", all(P$kA + P$kB == P$total), "\n")
cat("  frac_A == kA/total : ", all(abs(P$frac_A - P$kA/P$total) < 1e-9), "\n")
cat("  all p_adj < 0.05   : ", all(P$p_adj < 0.05), "\n")

P$dom_k <- pmax(P$kA, P$kB); P$sub_k <- pmin(P$kA, P$kB)
P$considered <- P$total + P$n_1to2
P$ret_dom <- (P$dom_k + P$n_1to2) / P$considered
P$ret_sub <- (P$sub_k + P$n_1to2) / P$considered
P$gap <- P$ret_dom - P$ret_sub
P$dom_chr <- ifelse(P$kA >= P$kB, P$chrA, P$chrB)
P$sub_chr <- ifelse(P$kA >= P$kB, P$chrB, P$chrA)

hr("1. RETENTION PER CHROMOSOME")
print(P[order(-P$gap), c("exp_pair","dom_chr","sub_chr","dom_k","sub_k","n_1to2",
                         "considered","ret_dom","ret_sub","gap")], row.names=FALSE)
cat(sprintf("\n  dominant members : %.3f - %.3f\n", min(P$ret_dom), max(P$ret_dom)))
cat(sprintf("  subordinate      : %.3f - %.3f\n", min(P$ret_sub), max(P$ret_sub)))
sep <- min(P$ret_dom) - max(P$ret_sub)
cat(sprintf("  separation between the two groups: %+.3f  %s\n", sep,
    if (sep > 0) "-- NO OVERLAP" else "-- THEY OVERLAP, the two-band claim fails"))
if (sep <= 0) cat("  !! the plot will not show clean bands; read the table above\n")
cat(sprintf("  informative genes %s | kept by both %s | considered %s\n",
    format(sum(P$total), big.mark=","), format(sum(P$n_1to2), big.mark=","),
    format(sum(P$considered), big.mark=",")))
write.csv(P[order(-P$gap),], file.path(OUT,"PFD_retention.csv"), row.names=FALSE)

hr("2. NULL: SAME COUNTS, LOSS AT RANDOM")
NS <- 4000
nul <- t(sapply(seq_len(nrow(P)), function(i){
  k  <- rbinom(NS, P$total[i], 0.5)
  hi <- pmax(k, P$total[i]-k); lo <- P$total[i]-hi
  g  <- (hi - lo)/P$considered[i]
  c(med_dom = median((hi + P$n_1to2[i])/P$considered[i]),
    med_sub = median((lo + P$n_1to2[i])/P$considered[i]),
    gap_med = median(g), gap_p95 = quantile(g, 0.95)) }))
for (i in seq_len(nrow(P)))
  cat(sprintf("    %-6s observed gap %.3f | null gap median %.3f  95th %.3f  %s\n",
      P$exp_pair[i], P$gap[i], nul[i,"gap_med"], nul[i,"gap_p95.95%"],
      if (P$gap[i] > nul[i,"gap_p95.95%"]) "beyond null" else "WITHIN NULL"))
cat(sprintf("  observed gap exceeds its own null 95th on %d of %d pairs\n",
    sum(P$gap > nul[,"gap_p95.95%"]), nrow(P)))
cat(sprintf("  observed gaps %.0f%%-%.0f%% | null gaps %.0f%%-%.0f%%\n",
    100*min(P$gap), 100*max(P$gap), 100*min(nul[,"gap_med"]), 100*max(nul[,"gap_med"])))

hr("3. PANELS")
p_votes <- function(){
  par(mar=c(0.5,0.5,3.4,0.5))
  plot(NA, xlim=c(0,100), ylim=c(0,100), axes=FALSE, xlab="", ylab="")
  mtext("Each ancestral gene votes once", side=3, adj=0, line=1.3, cex=1.05, font=2)
  text(0, 93, "is it still there on each copy of the chromosome?", adj=0, cex=0.84, col=GRY)
  rows <- list(list(y=74, a=TRUE,  b=TRUE,  t="kept on both",       s="no vote", c=GRY),
               list(y=56, a=TRUE,  b=FALSE, t="kept on one only",   s="a vote for that copy", c=WIN),
               list(y=38, a=FALSE, b=TRUE,  t="kept on the other",  s="a vote for that copy", c=LOSE),
               list(y=20, a=FALSE, b=FALSE, t="gone from both",     s="no vote", c=GRY))
  for (r in rows) {
    if (r$a) rect(2,r$y-5.5,17,r$y+5.5, col=WIN,  border=NA)
    else     rect(2,r$y-5.5,17,r$y+5.5, col=NA,   border=WIN,  lty=3)
    if (r$b) rect(20,r$y-5.5,35,r$y+5.5, col=LOSE, border=NA)
    else     rect(20,r$y-5.5,35,r$y+5.5, col=NA,   border=LOSE, lty=3)
    text(41, r$y+2.8, r$t, adj=0, cex=0.88)
    text(41, r$y-4.6, r$s, adj=0, cex=0.76, col=r$c) }
  text(2, 5, "solid = present     dashed = lost", adj=0, cex=0.74, col=GRY) }

p_bands <- function(){
  par(mar=c(4.4,5.4,3.6,6.2), xpd=NA)
  o  <- order(P$exp_pair); n <- nrow(P)
  yl <- c(min(c(P$ret_sub, nul[,"med_sub"]))-0.06, max(c(P$ret_dom, nul[,"med_dom"]))+0.05)
  plot(NA, xlim=c(0.4, n+0.6), ylim=yl, axes=FALSE, xlab="", ylab="")
  rect(0.4, min(P$ret_dom)-0.012, n+0.6, max(P$ret_dom)+0.012,
       col=adjustcolor(WIN, alpha.f=0.10), border=NA)
  rect(0.4, min(P$ret_sub)-0.012, n+0.6, max(P$ret_sub)+0.012,
       col=adjustcolor(LOSE, alpha.f=0.10), border=NA)
  for (j in seq_len(n)) { k <- o[j]
    segments(j+0.16, nul[k,"med_sub"], j+0.16, nul[k,"med_dom"], col=NUL, lwd=1.4)
    points(c(j+0.16, j+0.16), c(nul[k,"med_dom"], nul[k,"med_sub"]), pch=19, cex=0.9, col=NUL)
    segments(j-0.16, P$ret_sub[k], j-0.16, P$ret_dom[k], col="#B4B2A9", lwd=1.4)
    points(j-0.16, P$ret_dom[k], pch=19, cex=1.9, col=WIN)
    points(j-0.16, P$ret_sub[k], pch=19, cex=1.9, col=LOSE) }
  ax <- pretty(yl, 5); ax <- ax[ax >= yl[1] & ax <= yl[2]]
  axis(2, at=ax, labels=sprintf("%.0f%%", 100*ax), las=1, cex.axis=0.85, lwd=0.6)
  axis(1, at=seq_len(n), labels=sub("^chr","",P$exp_pair[o]), cex.axis=0.85, lwd=0.6)
  mtext("chromosome pair", side=1, line=2.5, cex=0.88)
  mtext("share of ancestral genes still present", side=2, line=3.8, cex=0.88)
  mtext("Sixteen chromosomes, two kinds", side=3, adj=0, line=2.0, cex=1.05, font=2)
  mtext(sprintf("no chromosome falls between the groups | %s ancestral genes",
        format(sum(P$considered), big.mark=",")), side=3, adj=0, line=0.7, cex=0.8, col=GRY)
  xr <- n+0.75
  text(xr, mean(range(P$ret_dom)), "kept\nmore", adj=0, cex=0.8, col="#0F6E56")
  text(xr, mean(range(P$ret_sub)), "kept\nfewer", adj=0, cex=0.8, col="#993C1D")
  text(xr, mean(c(nul[o[n],"med_dom"], nul[o[n],"med_sub"])), "if loss\nwere random",
       adj=0, cex=0.72, col=NUL) }

emit("PFD_1_votes", 5.2, 4.4, p_votes)
emit("PFD_2_bands", 6.6, 4.4, p_bands)
emit("PFD_ALL", 12.2, 4.6, function(){
  layout(matrix(1:2, nrow=1), widths=c(1, 1.32)); p_votes(); p_bands() })

hr("DONE")
print(list.files(OUT, pattern="^PFD"))
