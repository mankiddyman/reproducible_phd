#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Is the X/Y asymmetry gene loss or annotation quality?
# Anchor orthologs give a per-HOG expected protein length -> completeness ratio.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
res <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE)

part <- res %>% filter(anchor == "Nepenthes_gracilis") %>%
  transmute(pair_lab, X = winner, Y = ifelse(winner == chrA, chrB, chrA))
xy <- bind_rows(transmute(part, chr = X, side = "X", pair_lab),
                transmute(part, chr = Y, side = "Y", pair_lab))
cat("=== partition under test ===\n"); print(as.data.frame(part))

## ---- TEST 1: array-collapse rate ------------------------------------------
## A gene split into adjacent same-HOG models is called an array by GENESPACE
## and collapsed. Systematically more collapsing on Y = fragmentation.
arr <- bed %>% filter(genome == "Dionaea_muscipula") %>%
  group_by(chr) %>%
  summarise(n_models = n(), n_rep = sum(as.logical(isArrayRep), na.rm = TRUE),
            .groups = "drop") %>%
  mutate(frac_collapsed = (n_models - n_rep) / n_models) %>%
  inner_join(xy, by = "chr")
a2 <- arr %>% select(pair_lab, side, frac_collapsed) %>%
  pivot_wider(names_from = side, values_from = frac_collapsed)
cat("\n=== TEST 1: fraction of gene models collapsed as array members ===\n")
print(as.data.frame(a2), digits = 3)
cat(sprintf("Y > X in %d/8 | paired Wilcoxon p=%.3f  (Y higher = fragmentation on Y)\n",
            sum(a2$Y > a2$X), suppressWarnings(wilcox.test(a2$X, a2$Y, paired = TRUE)$p.value)))

## ---- completeness ratio ----------------------------------------------------
b <- bed %>% mutate(isRep = as.logical(isArrayRep)) %>% filter(is.na(isRep) | isRep)

ref <- b %>% filter(genome != "Dionaea_muscipula") %>%
  group_by(globHOG) %>%
  summarise(ref_len = median(pepLen, na.rm = TRUE), n_ref = n(), .groups = "drop") %>%
  filter(n_ref >= 2, ref_len > 0)

dio <- b %>% filter(genome == "Dionaea_muscipula") %>%
  inner_join(xy, by = "chr") %>% inner_join(ref, by = "globHOG") %>%
  mutate(lr = pepLen / ref_len)

st <- dio %>% group_by(globHOG) %>%
  summarise(n_dio = n(), n_pair = n_distinct(pair_lab), n_side = n_distinct(side),
            .groups = "drop")
dio <- left_join(dio, st, by = "globHOG")

cat(sprintf("\nDionaea genes with >=2 anchor orthologs: %d\n", nrow(dio)))

## ---- TEST 2: paired completeness, genes retained on BOTH homeologs ---------
pw <- dio %>% filter(n_dio == 2, n_pair == 1, n_side == 2) %>%
  select(globHOG, pair_lab, side, lr) %>%
  pivot_wider(names_from = side, values_from = lr) %>% filter(!is.na(X), !is.na(Y))
cat("\n=== TEST 2: completeness ratio, SAME gene on X vs Y (retained pairs) ===\n")
pp <- pw %>% group_by(pair_lab) %>%
  summarise(n = n(), medX = median(X), medY = median(Y),
            med_diff = median(X - Y),
            p = suppressWarnings(wilcox.test(X, Y, paired = TRUE)$p.value),
            .groups = "drop") %>% mutate(p_adj = p.adjust(p, "BH"))
print(as.data.frame(pp), digits = 3)
cat(sprintf("ALL pairs pooled: n=%d medX=%.3f medY=%.3f  paired Wilcoxon p=%.3g\n",
            nrow(pw), median(pw$X), median(pw$Y),
            suppressWarnings(wilcox.test(pw$X, pw$Y, paired = TRUE)$p.value)))

## ---- TEST 3: completeness of single-copy survivors, X-side vs Y-side ------
sc <- dio %>% filter(n_dio == 1)
cat("\n=== TEST 3: completeness ratio of single-copy survivors ===\n")
s3 <- sc %>% group_by(pair_lab, side) %>%
  summarise(n = n(), med_lr = median(lr), frac_short = mean(lr < 0.7), .groups = "drop") %>%
  pivot_wider(names_from = side, values_from = c(n, med_lr, frac_short))
print(as.data.frame(s3), digits = 3)
cat(sprintf("med_lr X>Y in %d/8 | Wilcoxon p=%.3f\n",
            sum(s3$med_lr_X > s3$med_lr_Y),
            suppressWarnings(wilcox.test(s3$med_lr_X, s3$med_lr_Y, paired = TRUE)$p.value)))
cat(sprintf("pooled single-copy: medX=%.3f (n=%d)  medY=%.3f (n=%d)  p=%.3g\n",
            median(sc$lr[sc$side=="X"]), sum(sc$side=="X"),
            median(sc$lr[sc$side=="Y"]), sum(sc$side=="Y"),
            suppressWarnings(wilcox.test(lr ~ side, data = sc)$p.value)))

## ---- TEST 4: retention re-scored on complete models only ------------------
## If Y's "loss" is fragmentation, dropping incomplete models should NOT move
## frac_A; if fragmentation inflates Y's gene count with junk, it will.
rescore <- function(thr) {
  keep <- dio %>% filter(n_dio == 1, lr >= thr)
  part %>% rowwise() %>%
    mutate(kX = sum(keep$side == "X" & keep$pair_lab == pair_lab),
           kY = sum(keep$side == "Y" & keep$pair_lab == pair_lab)) %>%
    ungroup() %>% mutate(frac_X = kX / (kX + kY), thr = thr) %>%
    select(pair_lab, thr, kX, kY, frac_X)
}
rs <- bind_rows(lapply(c(0, 0.5, 0.7, 0.9), rescore)) %>%
  pivot_wider(names_from = thr, values_from = c(kX, kY, frac_X))
cat("\n=== TEST 4: retention after dropping incomplete models ===\n")
cat("frac_X = X's share of single-copy HOGs. Stable across thresholds = robust.\n")
print(as.data.frame(select(rs, pair_lab, starts_with("frac_X"))), digits = 3)

write_csv(pw, "completeness_paired.csv")
write_csv(rs, "retention_by_completeness.csv")
