#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ape); library(ragg) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

meta <- read_tsv("full14/tip_meta.tsv", show_col_types = FALSE)
fullset <- readLines("full14/list_full.txt")
meta <- meta %>% mutate(
  lab = ifelse(genome == "Nepenthes_gracilis", "NEPENTHES",
        paste0(ifelse(grepl("Dionaea", genome), "Dio",
               ifelse(grepl("regia", genome), "reg", "bin")), " ",
               sub("_hap1$|_collapsed$", "", sub("^chr", "chr", chr)))))
key <- function(x) gsub("@", "_", gsub("['\"]", "", x))

draw <- function(a, cex_tip = 1.15) {
  f <- file.path("full14/tre", paste0(a, ".treefile"))
  if (!file.exists(f)) { plot.new(); title(a); return(invisible()) }
  tr <- read.tree(f)
  m <- meta[match(key(tr$tip.label), key(meta$tip)), ]
  nep <- tr$tip.label[which(m$genome == "Nepenthes_gracilis")]
  if (length(nep) == 1) tr <- root(tr, nep, resolve.root = TRUE)
  m <- meta[match(key(tr$tip.label), key(meta$tip)), ]
  # IQ-TREE -B writes "SHaLRT/UFboot" or just UFboot; keep UFboot
  sup <- tr$node.label
  if (!is.null(sup)) {
    sup <- ifelse(grepl("/", sup), sub(".*/", "", sup), sup)
    sup <- suppressWarnings(as.numeric(sup))
    sup[is.na(sup) | sup < 50] <- NA
  }
  tr$tip.label <- m$lab
  col <- ifelse(grepl("^Dio", tr$tip.label), "#C0392B",
         ifelse(grepl("^reg", tr$tip.label), "#5B4EA8",
         ifelse(grepl("^bin", tr$tip.label), "#12795E", "#333333")))
  fnt <- ifelse(grepl("^Dio", tr$tip.label), 2, 1)
  plot(tr, tip.color = col, font = fnt, cex = cex_tip, edge.width = 2.2,
       label.offset = 0.012, no.margin = FALSE, x.lim = NULL)
  if (!is.null(sup) && any(!is.na(sup)))
    nodelabels(ifelse(is.na(sup), "", sup), frame = "none", cex = cex_tip * 0.72,
               adj = c(1.25, -0.55), col = "#555555")
  add.scale.bar(lwd = 1.6, cex = cex_tip * 0.7)
  title(sub("\\.t1$", "", a), cex.main = cex_tip * 0.95, font.main = 1)
}

## one page, 14 panels, large
agg_png("FIG7_full14_trees.png", width = 4200, height = 3000, res = 260)
par(mfrow = c(4, 4), mar = c(2.2, 0.6, 2.6, 0.6), xpd = NA)
for (a in fullset) draw(a, cex_tip = 1.05)
plot.new(); legend("center", bty = "n", cex = 1.5, text.font = 2,
  legend = c("Dionaea", "D. regia", "D. binata", "Nepenthes"),
  text.col = c("#C0392B", "#5B4EA8", "#12795E", "#333333"))
dev.off()

## and one big panel per tree, for reading detail
dir.create("full14/png", showWarnings = FALSE)
for (a in fullset) {
  agg_png(file.path("full14/png", paste0(a, ".png")), width = 2000, height = 1400, res = 220)
  par(mar = c(3, 1, 3, 1), xpd = NA)
  draw(a, cex_tip = 1.5)
  dev.off()
}
cat(sprintf("WROTE FIG7_full14_trees.png (4200x3000) and %d single-tree PNGs in full14/png/\n",
            length(fullset)))
