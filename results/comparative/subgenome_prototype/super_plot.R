#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
suppressPackageStartupMessages({ library(ape); library(ragg); library(dplyr) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

fs <- list.files("super", pattern = "^chr[0-9]+_dom\\.treefile$", full.names = TRUE)
cat(sprintf("region trees: %d\n", length(fs)))

for (f in fs) {
  r <- sub("\\.treefile$", "", basename(f))
  tr <- read.tree(f)
  if ("NEP" %in% tr$tip.label) tr <- root(tr, "NEP", resolve.root = TRUE)
  col <- ifelse(grepl("^Dio", tr$tip.label), "#C0392B",
         ifelse(grepl("^reg", tr$tip.label), "#5B4EA8",
         ifelse(grepl("^bin", tr$tip.label), "#12795E", "#111111")))
  h <- max(1400, 45 * Ntip(tr))
  agg_png(sprintf("FIG8_%s.png", r), width = 2400, height = h, res = 200)
  par(mar = c(3, 1, 3, 1), xpd = NA)
  plot(tr, tip.color = col, cex = 1.0, edge.width = 2, label.offset = 0.004,
       font = ifelse(grepl("^Dio", tr$tip.label), 2, 1))
  if (!is.null(tr$node.label)) {
    s <- suppressWarnings(as.numeric(sub(".*/", "", tr$node.label)))
    nodelabels(ifelse(is.na(s) | s < 70, "", s), frame = "none",
               cex = 0.75, adj = c(1.2, -0.5), col = "#666666")
  }
  add.scale.bar(lwd = 2)
  title(sprintf("%s  (%d blocks)", r, Ntip(tr)), cex.main = 1.2, font.main = 1)
  dev.off()

  sf <- sprintf("super/%s_scf.cf.stat", r)
  if (file.exists(sf)) {
    d <- read.table(sf, header = TRUE, comment.char = "#")
    cat(sprintf("%-12s tips=%-3d median sCF=%.1f  branches sCF>50: %d/%d\n",
                r, Ntip(tr), median(d$sCF, na.rm = TRUE),
                sum(d$sCF > 50, na.rm = TRUE), nrow(d)))
  }
}
cat("\nWROTE: FIG8_chr*_dom.png\n")
