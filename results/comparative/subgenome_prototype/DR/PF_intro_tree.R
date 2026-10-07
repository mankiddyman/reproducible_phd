#!/usr/bin/env Rscript
# ============================================================================
# PF_intro_tree -- the study system, for the introduction
#
# Seven tips: Nepenthes, Dionaea, and the five Drosera. Topology from
# DR/tree/astral/astral_all.tre collapsed to species (both subgenome tips of a
# species map to one tip). Thick branch = holocentric, thin = monocentric.
# Core Drosera bracketed; regia sits outside it, with Dionaea.
#
# Centromere assignment is from unpublished lab work -- change CENT below if it
# changes. No dates, no support values: this is orientation, not a result.
#
# OUT presentation_figures/PFintro_tree.{pdf,png}
# ============================================================================
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=200)
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

HOLO <- c("regia","scorpioides","paradoxa")
MONO <- c("capensis","binata")
LW_HOLO <- 4.6; LW_MONO <- 1.6
HOL_C <- "#534AB7"; MON_C <- "#888780"; N_COL <- "#6E6D69"; BR <- "#BA7517"

# tips top to bottom, with the x at which each terminal branch starts
TIP <- list(
  list(y=7, x=30, g="Nepenthes",   it="Nepenthes gracilis",  cent=NA),
  list(y=6, x=46, g="Dionaea",     it="Dionaea muscipula",   cent=NA),
  list(y=5, x=46, g="regia",       it="Drosera regia",       cent="holo"),
  list(y=4, x=46, g="capensis",    it="Drosera capensis",    cent="mono"),
  list(y=3, x=58, g="binata",      it="Drosera binata",      cent="mono"),
  list(y=2, x=70, g="paradoxa",    it="Drosera paradoxa",    cent="holo"),
  list(y=1, x=70, g="scorpioides", it="Drosera scorpioides", cent="holo"))
TX <- 82

draw <- function(){
  par(mar=c(0.6,0.6,3.0,0.6))
  plot(NA, xlim=c(0,170), ylim=c(0.1,8.1), axes=FALSE, xlab="", ylab="")
  sg <- function(x1,y1,x2,y2,lwd=1.6,col="black")
    segments(x1,y1,x2,y2,lwd=lwd,lend=1,col=col)

  sg(2,4.75,10,4.75)                      # root stub
  sg(10,2.5,10,7)                         # Nepenthes vs the rest
  sg(10,7,30,7)
  sg(10,2.5,30,2.5)                       # crown Droseraceae
  sg(30,1.75,30,5.5)
  sg(30,5.5,46,5.5)                       # (Dionaea, regia)
  sg(46,5,46,6)
  sg(30,1.75,46,1.75)                     # core Drosera
  sg(46,1.5,46,4)
  sg(46,1.5,58,1.5)
  sg(58,1.5,58,3)
  sg(58,1.5,70,1.5)
  sg(70,1,70,2)

  for (t in TIP) {
    lw <- if (identical(t$cent,"holo")) LW_HOLO else LW_MONO
    cl <- if (is.na(t$cent)) "black" else if (t$cent=="holo") HOL_C else MON_C
    sg(t$x, t$y, TX, t$y, lwd=lw, col=cl)
    text(TX+3, t$y, t$it, adj=0, cex=1.15, font=3) }

  ytop <- 4; ybot <- 1
  sg(160, ybot, 160, ytop, lwd=1.2, col=BR)
  sg(157, ybot, 160, ybot, lwd=1.2, col=BR)
  sg(157, ytop, 160, ytop, lwd=1.2, col=BR)
  text(164, (ytop+ybot)/2, "core", adj=0, cex=1.15, font=2, col=BR)
  text(164, (ytop+ybot)/2 - 0.42, "Drosera", adj=0, cex=1.15, font=2, col=BR)

  sg(10, 0.35, 22, 0.35, lwd=LW_HOLO, col=HOL_C); text(24, 0.35, "holocentric", adj=0, cex=1.0)
  sg(62, 0.35, 74, 0.35, lwd=LW_MONO, col=MON_C); text(76, 0.35, "monocentric", adj=0, cex=1.0)

  mtext("The system", side=3, adj=0, line=1.3, cex=1.3, font=2)
  mtext("centromere type varies within the genus", side=3, adj=0, line=0.1,
        cex=0.92, col=N_COL) }

hr("RENDER")
cat("  holocentric:", paste(HOLO, collapse=", "), "\n")
cat("  monocentric:", paste(MONO, collapse=", "), "\n")
W <- 9.2; H <- 5.2
pdf(file.path(OUT,"PFintro_tree.pdf"), width=W, height=H, useDingbats=FALSE); draw(); dev.off()
png(file.path(OUT,"PFintro_tree.png"), width=W*160, height=H*160, res=160); draw(); dev.off()
cat("  wrote PFintro_tree.pdf / .png\n")
