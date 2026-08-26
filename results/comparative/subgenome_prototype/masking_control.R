#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Does soft-masked fraction explain the X/Y gene-count gap?
# If genes-per-unmasked-Mb equalises, the "fractionation" is a masking artifact.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

msk <- read.table("masking_by_chr.tsv", col.names = c("chr","lc","tot"),
                  stringsAsFactors = FALSE) %>% mutate(frac_masked = lc/tot)
bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
res <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE)

part <- res %>% filter(anchor == "Nepenthes_gracilis") %>%
  transmute(pair_lab, X = winner, Y = ifelse(winner == chrA, chrB, chrA),
            frac_X_obs = ifelse(winner == chrA, frac_A, 1 - frac_A))

gc <- bed %>% filter(genome == "Dionaea_muscipula") %>%
  mutate(isRep = as.logical(isArrayRep)) %>% filter(is.na(isRep) | isRep) %>%
  count(chr, name = "n_genes") %>% inner_join(msk, by = "chr") %>%
  mutate(unmasked_Mb = (tot - lc)/1e6, total_Mb = tot/1e6,
         dens_total = n_genes/total_Mb, dens_unmasked = n_genes/unmasked_Mb)

cat("=== per-chromosome masking and density ===\n")
print(as.data.frame(gc %>% select(chr, total_Mb, frac_masked, n_genes,
                                  dens_total, dens_unmasked)), digits = 3)

cmp <- part %>%
  left_join(select(gc, chr, mX = frac_masked, gX = n_genes, uX = unmasked_Mb,
                   dtX = dens_total, duX = dens_unmasked), by = c("X"="chr")) %>%
  left_join(select(gc, chr, mY = frac_masked, gY = n_genes, uY = unmasked_Mb,
                   dtY = dens_total, duY = dens_unmasked), by = c("Y"="chr")) %>%
  mutate(gene_ratio = gX/gY,
         exp_from_unmasked = uX/(uX+uY),
         obs_gene_share = gX/(gX+gY))

cat("\n=== X vs Y ===\n")
print(as.data.frame(select(cmp, pair_lab, mX, mY, gX, gY, duX, duY,
                           exp_from_unmasked, obs_gene_share, frac_X_obs)), digits = 3)

cat(sprintf("\nmasked  Y>X in %d/8 | paired Wilcoxon p=%.3f\n",
            sum(cmp$mY > cmp$mX), suppressWarnings(wilcox.test(cmp$mX, cmp$mY, paired=TRUE)$p.value)))
cat(sprintf("dens_total     X>Y in %d/8 | p=%.3f\n", sum(cmp$dtX > cmp$dtY),
            suppressWarnings(wilcox.test(cmp$dtX, cmp$dtY, paired=TRUE)$p.value)))
cat(sprintf("dens_UNMASKED  X>Y in %d/8 | p=%.3f   <- if this collapses, masking explains it\n",
            sum(cmp$duX > cmp$duY),
            suppressWarnings(wilcox.test(cmp$duX, cmp$duY, paired=TRUE)$p.value)))
ct <- cor.test(cmp$exp_from_unmasked, cmp$obs_gene_share, method = "spearman")
cat(sprintf("Spearman(unmasked-bp share, gene share) rho=%.3f p=%.4f\n", ct$estimate, ct$p.value))
write_csv(cmp, "masking_control.csv")
