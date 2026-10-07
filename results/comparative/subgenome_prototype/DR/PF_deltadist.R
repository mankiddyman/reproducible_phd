#!/usr/bin/env Rscript
# ============================================================================
# PF_deltadist -- per-gene Delta by species, ordered by distance from Dionaea
#
# Panels run regia -> core Drosera. Delta for an A-derived copy is
# T_spec/T_AB - 1, so a lineage that split from Dionaea more recently has a
# smaller T_spec and a Delta further from zero: the separation SHARPENS with
# proximity to Dionaea. regia being sister to Dionaea therefore predicts its
# bimodality, independently of the tree.
#
# The A-side fraction is not 67% because per-gene misassignment e pulls it in:
#   observed = 0.67(1-e) + 0.33e = 0.67 - 0.34e
# so e is recovered per species as (0.67 - observed)/0.34 and printed.
#
# Centromere type is a marker beside each species name, not the grouping --
# it cuts across the phylogeny.
#
# IN   DR/out/DR02_gene_delta.csv
# OUT  presentation_figures/PFdelta_bySpecies.{pdf,png}
#      presentation_figures/PFdelta_split.csv
# ============================================================================
suppressMessages(library(data.table))
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }

A_COL <- "#0F6E56"; B_COL <- "#993C1D"; N_COL <- "#6E6D69"
A_FIL <- "#1D9E75"; B_FIL <- "#D85A30"

HOLO  <- c("Drosera_regia","Drosera_scorpioides","Drosera_paradoxa")
ORDER <- c("Drosera_regia","Drosera_capensis","Drosera_scorpioides",
           "Drosera_binata","Drosera_paradoxa")     # regia first, then core
BW    <- 0.10

hr("0. LOAD")
GD <- "DR/out/DR02_gene_delta.csv"
if (!file.exists(GD)) die("%s not found", GD)
G <- fread(GD)[is.finite(delta)]
sp <- unique(G$genome)
cat("  copies:", format(nrow(G), big.mark=","), "| species:", paste(sp, collapse=", "), "\n")
miss <- setdiff(sp, ORDER); if (length(miss)) die("not in ORDER: %s", paste(miss, collapse=", "))
ORDER <- ORDER[ORDER %in% sp]
G[, cent := fifelse(genome %in% HOLO, "holocentric", "monocentric")]
cat("  holocentric:", paste(sub("Drosera_","D. ", HOLO), collapse=", "), "\n")

hr("1. SPLIT, AND THE IMPLIED PER-GENE ERROR")
S <- G[, .(copies=.N, frac_A=mean(delta<0), median=median(delta), sd=sd(delta),
           strong=mean(abs(delta)>0.5)), by=genome]
S[, err_implied := (2/3 - frac_A)/(1/3)]
S[, ratio := sprintf("%.2f : 1", frac_A/(1-frac_A))]
S <- S[match(ORDER, genome)]
print(S, row.names=FALSE)
cat(sprintf("\n  pooled frac_A %.3f (%.2f : 1) | AAB predicts 0.667\n",
    mean(G$delta<0), mean(G$delta<0)/(1-mean(G$delta<0))))
cat("  a true 2:1 seen through per-gene error e gives 0.667 - 0.34e;\n")
cat("  e recovered per species above. DR02's noise model predicts e = 0.31.\n")
fwrite(S, file.path(OUT,"PFdelta_split.csv"))

hr("2. FIGURE")
XL <- c(-1,1)
dens <- function(v) density(v, bw=BW, from=XL[1], to=XL[2], n=512)

panel <- function(g){
  v <- G[genome==g, delta]; d <- dens(v); ym <- max(d$y)
  plot(NA, xlim=XL, ylim=c(0, ym*1.30), axes=FALSE, xlab="", ylab="")
  li <- d$x<=0; ri <- d$x>=0
  polygon(c(d$x[li],0), c(d$y[li],0), col=adjustcolor(A_FIL, alpha.f=0.30), border=NA)
  polygon(c(0,d$x[ri]), c(0,d$y[ri]), col=adjustcolor(B_FIL, alpha.f=0.30), border=NA)
  lines(d$x, d$y, lwd=2.0, col="black")
  segments(0, 0, 0, ym*1.06, lwd=1.1, col="black")
  axis(1, at=c(-1,-0.5,0,0.5,1), labels=c("\u22121","","0","","+1"),
       cex.axis=1.0, lwd=0.6, padj=-0.3)
  f <- mean(v<0)
  text(-0.95, ym*1.20, sprintf("%.0f%%", 100*f),     adj=0, cex=1.35, font=2, col=A_COL)
  text( 0.95, ym*1.20, sprintf("%.0f%%", 100*(1-f)), adj=1, cex=1.35, font=2, col=B_COL)
  hol <- G[genome==g, cent][1] == "holocentric"
  points(-0.95, ym*1.02, pch=if (hol) 19 else 1, cex=1.0, col=N_COL)
  text(-0.88, ym*1.02, sub("Drosera_","D. ", g), adj=0, cex=1.15, font=3) }

inset_tree <- function(){
  plot(NA, xlim=c(0,100), ylim=c(100,0), axes=FALSE, xlab="", ylab="")
  segments(c(6,16,16,38,38,16,52,52), c(50,26,26,10,10,74,74,58),
           c(16,16,38,52,52,52,52,74), c(50,74,26,10,42,74,58,58),
           lwd=1.8, col="black")
  text(56, 14, "Dionaea",  adj=0, cex=1.0, font=3, col=N_COL)
  text(56, 46, "D. regia", adj=0, cex=1.0, font=3)
  text(78, 62, "core",     adj=0, cex=1.0, font=3)
  text(78, 76, "Drosera",  adj=0, cex=1.0, font=3)
  text(0, 96, "regia is sister to Dionaea, not to core Drosera",
       adj=0, cex=0.92, col=N_COL) }

emit <- function(nm,w,h,fn){
  pdf(file.path(OUT,paste0(nm,".pdf")), width=w, height=h, useDingbats=FALSE); fn(); dev.off()
  png(file.path(OUT,paste0(nm,".png")), width=w*160, height=h*160, res=160); fn(); dev.off()
  cat("  wrote ", nm, "\n", sep="") }

emit("PFdelta_bySpecies", 12, 5.6, function(){
  par(mfrow=c(2,3), mar=c(3.4,1.0,1.4,1.0), oma=c(4.0,1.0,4.2,1.0))
  for (g in ORDER) panel(g)
  if (length(ORDER) < 6) { par(mar=c(3.4,2.0,1.4,1.0)); inset_tree() }
  mtext("\u0394    (\u22121 = pure A,  +1 = pure B)", side=1, outer=TRUE, line=2.0, cex=1.05)
  mtext("Every species splits the same way \u2014 and the closer to Dionaea, the sharper",
        side=3, outer=TRUE, adj=0, line=2.2, cex=1.2, font=2)
  mtext(sprintf("%s gene copies | filled dot = holocentric, open = monocentric | a true 2:1 seen through per-gene noise reads as ~56\u201362%%",
        format(nrow(G), big.mark=",")),
        side=3, outer=TRUE, adj=0, line=0.6, cex=0.88, col=N_COL) })

hr("DONE")
print(list.files(OUT, pattern="^PFdelta_bySpecies"))
