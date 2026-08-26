#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Block-level topology. NO use of the X/Y partition anywhere -- Dionaea has HE too,
# so its chromosome labels are locally wrong. Blocks are the unit for every species.
# Per ancestral region: aggregate gene-tree topology into a block x block affinity
# matrix (MRCA clade size, rank-normalised per tree -> rate-free), build a block tree,
# then compare shapes across the 8 regions.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ape); library(ggplot2); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
MINGENES <- 5

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
blk <- read_csv(file.path(GSD, "syntenicBlock_coordinates.csv"), show_col_types = FALSE)
meta <- read_tsv("wgd/tip_meta.tsv", show_col_types = FALSE)
anch <- read_csv("synteny_ortho_hits.csv", show_col_types = FALSE) %>% distinct(nep_gene, nep_chr)

nb <- blk %>% filter(genome1 == "Nepenthes_gracilis" | genome2 == "Nepenthes_gracilis") %>%
  mutate(n1 = genome1 == "Nepenthes_gracilis",
         sp = ifelse(n1, genome2, genome1), nep_chr = ifelse(n1, chr1, chr2),
         sp_chr = ifelse(n1, chr2, chr1),
         sp_s = ifelse(n1, pmin(startBp2, endBp2), pmin(startBp1, endBp1)),
         sp_e = ifelse(n1, pmax(startBp2, endBp2), pmax(startBp1, endBp1)),
         nhit = ifelse(n1, nHits2, nHits1)) %>%
  filter(sp != "Nepenthes_gracilis", grepl("_dom$", nep_chr), !is.na(sp_s)) %>%
  mutate(block = paste0(sub("^Dros?era_|^Dionaea_", "", sp), ":",
                        sub("_hap1$|_collapsed$|_sg[0-9]+_s[0-9]+$", "",
                            sub("^chr", "c", sp_chr)),
                        ifelse(grepl("_sg", sp_chr),
                               sub("^.*_(s[0-9]+)$", "-\\1", sp_chr), ""), "#", blkID)) %>%
  select(block, sp, nep_chr, sp_chr, sp_s, sp_e, nhit)
cat(sprintf("Nepenthes-anchored blocks: %d\n", nrow(nb)))

# FIX: tip_meta already has chr; do not re-join it
pos <- bed %>% transmute(gene = id, genome, gmid = (start + end) / 2, gord = ord)
tips <- meta %>% left_join(pos, by = c("gene", "genome")) %>%
  left_join(anch, by = "nep_gene")
stopifnot("chr" %in% names(tips), "gmid" %in% names(tips))

ab <- function(sp_, chr_, mid_, nchr_) {
  if (is.na(mid_) || is.na(nchr_)) return(NA_character_)
  c <- nb[nb$sp == sp_ & nb$sp_chr == chr_ & nb$nep_chr == nchr_ &
          nb$sp_s <= mid_ & nb$sp_e >= mid_, ]
  if (!nrow(c)) return(NA_character_)
  c$block[which.max(c$nhit)]
}
i <- which(tips$genome != "Nepenthes_gracilis")
tips$block <- NA_character_
tips$block[i] <- mapply(ab, tips$genome[i], tips$chr[i], tips$gmid[i], tips$nep_chr[i])
cat(sprintf("tips placed in a block: %.1f%% (%d/%d)\n",
            100*mean(!is.na(tips$block[i])), sum(!is.na(tips$block[i])), length(i)))
cat(sprintf("blocks with >=%d genes: %d\n", MINGENES,
            sum(table(tips$block) >= MINGENES)))

## ---- per gene tree: rank-normalised MRCA clade size for every tip pair -------
fs <- list.files("wgd/tre", pattern = "\\.tre$", full.names = TRUE)
one <- function(f) {
  tr <- tryCatch(read.tree(f), error = function(e) NULL)
  if (is.null(tr) || Ntip(tr) < 4) return(NULL)
  m <- tips[match(tr$tip.label, tips$tip), ]
  if (any(is.na(m$genome))) return(NULL)
  nep <- tr$tip.label[m$genome == "Nepenthes_gracilis"]
  if (length(nep) != 1) return(NULL)
  tr <- tryCatch(root(tr, nep, resolve.root = TRUE), error = function(e) NULL)
  if (is.null(tr)) return(NULL)
  m <- tips[match(tr$tip.label, tips$tip), ]
  k <- which(!is.na(m$block))
  if (length(k) < 2) return(NULL)
  nd <- node.depth(tr)
  cb <- combn(k, 2)
  d <- vapply(seq_len(ncol(cb)), function(z)
    nd[getMRCA(tr, tr$tip.label[cb[, z]])], 1)
  data.frame(region = m$nep_chr[k[1]],
             b1 = m$block[cb[1, ]], b2 = m$block[cb[2, ]],
             g1 = m$gene[cb[1, ]], g2 = m$gene[cb[2, ]],
             d = rank(d) / length(d),          # rank-normalised -> rate-free
             stringsAsFactors = FALSE)
}
v <- bind_rows(lapply(fs, one)) %>% filter(b1 != b2)
write_csv(v, "block_pair_affinity.csv")
cat(sprintf("block-pair observations: %d\n", nrow(v)))

## ---- within-block consistency: do genes in a block agree on nearest partner? --
near <- bind_rows(transmute(v, region, from = b1, to = b2, d),
                  transmute(v, region, from = b2, to = b1, d)) %>%
  mutate(sp_to = sub(":.*", "", to))
cons <- near %>% group_by(region, from, sp_to, to) %>%
  summarise(n = n(), md = mean(d), .groups = "drop") %>%
  group_by(region, from, sp_to) %>%
  summarise(n_obs = sum(n), best = to[which.min(md)],
            top_share = max(n) / sum(n),
            margin = round(diff(sort(md)[1:2]), 3), .groups = "drop") %>%
  filter(n_obs >= MINGENES)
write_csv(cons, "block_consistency.csv")
cat("\n=== within-block consistency (share of genes agreeing on modal partner) ===\n")
print(as.data.frame(cons %>% group_by(sp_from = sub(":.*", "", from), sp_to) %>%
  summarise(blocks = n(), median_share = round(median(top_share), 3), .groups = "drop")))

## ---- per-region block tree ---------------------------------------------------
agg <- v %>% group_by(region, b1, b2) %>%
  summarise(n = n(), d = mean(d), .groups = "drop") %>% filter(n >= 3)

regs <- sort(unique(agg$region))
agg_png("FIG6_block_trees.png", width = 1600, height = 950, res = 105)
par(mfrow = c(2, 4), mar = c(1, 1, 2.5, 1))
shapes <- list()
for (r in regs) {
  s <- filter(agg, region == r)
  bl <- sort(unique(c(s$b1, s$b2)))
  M <- matrix(NA_real_, length(bl), length(bl), dimnames = list(bl, bl))
  M[cbind(s$b1, s$b2)] <- s$d; M[cbind(s$b2, s$b1)] <- s$d; diag(M) <- 0
  keep <- bl[colSums(is.na(M)) == 0]
  if (length(keep) < 4) { plot.new(); title(paste0(r, "\n(too sparse)")); next }
  M <- M[keep, keep]
  tr <- tryCatch(nj(as.dist(M)), error = function(e) NULL)
  if (is.null(tr)) { plot.new(); title(r); next }
  tr <- midpoint_root <- tryCatch(root(tr, which.max(rowSums(M)), resolve.root = TRUE),
                                  error = function(e) tr)
  col <- ifelse(grepl("^muscipula", tr$tip.label), "#d95f02",
         ifelse(grepl("^regia", tr$tip.label), "#7570b3",
         ifelse(grepl("^binata", tr$tip.label), "#1b9e77", "grey50")))
  plot(tr, tip.color = col, cex = 0.7, edge.width = 1.3)
  title(sprintf("%s  (%d blocks)", r, length(keep)), cex.main = 1)
  # composition of the two halves of the deepest split
  ch <- tryCatch({
    d <- dist.nodes(tr); rt <- Ntip(tr) + 1
    kids <- tr$edge[tr$edge[, 1] == rt, 2]
    lapply(kids, function(k) if (k <= Ntip(tr)) tr$tip.label[k]
                             else extract.clade(tr, k)$tip.label)
  }, error = function(e) NULL)
  if (!is.null(ch) && length(ch) == 2)
    shapes[[r]] <- sapply(ch, function(x) paste(sort(table(sub(":.*", "", x))), collapse = "/"))
}
dev.off()

cat("\n=== composition of the two halves of each region's deepest split ===\n")
for (r in names(shapes)) cat(sprintf("%-10s  %s  ||  %s\n", r, shapes[[r]][1], shapes[[r]][2]))
cat("\nsame composition pattern across regions => consistent subgenome structure;\n")
cat("the ASYMMETRY (e.g. 1+1 vs 2+2) is what orients regions without any X/Y label\n")

## ---- HE detector: gene-by-gene partner flips along a block -------------------
big <- cons %>% filter(sp_to != sub(":.*", "", from)) %>%
  group_by(from) %>% slice_max(n_obs, n = 1) %>% ungroup() %>%
  slice_max(n_obs, n = 12)
det <- near %>% semi_join(select(big, from, sp_to), by = c("from", "sp_to")) %>%
  left_join(select(tips, gene, gord), by = c("from" = "gene"))
cat(sprintf("\nWROTE: FIG6_block_trees.png block_pair_affinity.csv block_consistency.csv\n"))
