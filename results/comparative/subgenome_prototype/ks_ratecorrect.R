#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Raw dS is rate x time, and Droseraceae rates differ between lineages, so raw modes are
# not comparable. Nepenthes split from all six at the same time, so dS(species->Nepenthes)
# measures that species' RATE against a constant. The ratio cancels the rate and gives a
# relative age: (time since the species' own copies split) / (time since Nepenthes split).
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ragg); library(patchwork)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
SPORD <- c("Dionaea_muscipula","Drosera_regia","Drosera_binata",
           "Drosera_paradoxa","Drosera_scorpioides","Drosera_capensis")
NEP <- "Nepenthes_gracilis"

k <- read_csv("ks/pairwise_ks.csv", show_col_types = FALSE) %>%
  filter(!is.na(dS), dS > 0, dS < 5, codons >= 100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types = FALSE)
cn <- tm %>% count(nep_gene, genome, name = "copies")

## ---- 1. the rate calibration: each species vs Nepenthes -------------------
NEPD <- bind_rows(
    k %>% filter(sp1 == NEP, sp2 != NEP) %>% transmute(anchor, species = sp2, dS),
    k %>% filter(sp2 == NEP, sp1 != NEP) %>% transmute(anchor, species = sp1, dS)) %>%
  group_by(anchor, species) %>% summarise(dS_nep = median(dS), .groups = "drop")

cat("=== RATE CALIBRATION: dS to Nepenthes (same true divergence time for all) ===\n")
rate <- NEPD %>% group_by(species) %>%
  summarise(n = n(), median_dS_nep = round(median(dS_nep), 3), .groups = "drop") %>%
  arrange(match(species, SPORD)) %>%
  mutate(rel_rate = round(median_dS_nep / median_dS_nep[species == "Dionaea_muscipula"], 2))
print(as.data.frame(rate), row.names = FALSE)
cat("rel_rate = how fast that lineage evolves relative to Dionaea.\n")
cat("If these differ, raw within-species dS CANNOT be compared across species.\n")

## ---- 2. rate-corrected within-species divergence --------------------------
W <- k %>% filter(sp1 == sp2, sp1 != NEP) %>% rename(species = sp1) %>%
  select(anchor, species, seq1, seq2, dS) %>%
  inner_join(NEPD, by = c("anchor", "species")) %>%
  left_join(cn, by = c("anchor" = "nep_gene", "species" = "genome")) %>%
  filter(!is.na(copies), dS_nep > 0.05) %>%
  mutate(ratio = dS / dS_nep,
         kgrp = ifelse(copies >= 4, "4+ copies", paste0(copies, " copies")),
         species = factor(species, SPORD))
cat(sprintf("\nwithin-species pairs with a Nepenthes calibration: %d\n", nrow(W)))

md <- function(x, to = 3) { if (length(x) < 20) return(NA_real_)
  d <- density(x, from = 0, to = to, n = 2048); d$x[which.max(d$y)] }

cat("\n=== TEST 1, rate-corrected: ratio = dS(own copies) / dS(to Nepenthes) ===\n")
res <- W %>% group_by(species, kgrp) %>%
  summarise(n = n(),
            raw_mode = round(md(dS, 2), 3),
            ratio_mode = round(md(ratio), 3),
            ratio_median = round(median(ratio), 3), .groups = "drop") %>%
  filter(n >= 20)
print(as.data.frame(res), row.names = FALSE)
cat("\nSAME ratio across species  => one shared progenitor pair; raw differences were rate\n")
cat("DIFFERENT ratios           => genuinely different divergence times, different progenitors\n")

d2 <- filter(W, copies == 2)
cat("\n=== 2-copy loci only (one comparison per locus, cleanest) ===\n")
print(as.data.frame(d2 %>% group_by(species) %>%
  summarise(n = n(), raw = round(md(dS,2),3), ratio = round(md(ratio),3),
            ratio_IQR = paste0(round(quantile(ratio,.25),2), "-", round(quantile(ratio,.75),2)),
            .groups="drop")), row.names = FALSE)

## ---- 3. what the near-zero spikes are --------------------------------------
cat("\n=== near-zero dS: recent duplicates or gene conversion, NOT homeologs ===\n")
print(as.data.frame(W %>% group_by(species) %>%
  summarise(n = n(), frac_dS_lt_0.1 = round(mean(dS < 0.1), 3),
            frac_saturated = round(mean(dS > 1), 2), .groups = "drop")), row.names = FALSE)

## ---- plots -----------------------------------------------------------------
ref <- md(filter(W, species == "Dionaea_muscipula", copies == 2)$ratio)
p1 <- ggplot(W, aes(ratio)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = .05, boundary = 0,
                 fill = "#12795E", colour = NA, alpha = .55) +
  geom_density(colour = "#0d5c47", linewidth = .7, bw = .06) +
  geom_vline(xintercept = ref, colour = "#C0392B", linetype = 2, linewidth = .6) +
  facet_grid(species ~ kgrp, scales = "free_y", switch = "y") +
  coord_cartesian(xlim = c(0, 2)) +
  labs(title = "TEST 1, RATE-CORRECTED — relative age of a species' own gene copies",
       subtitle = paste0("ratio = dS(own copies) / dS(that species to Nepenthes), which cancels lineage rate differences\n",
                         "dashed red = the Dionaea 2-copy value (", round(ref,3),
                         "). Same position across species = shared progenitor pair."),
       x = "dS(own copies) / dS(to Nepenthes)", y = NULL) +
  theme_bw(base_size = 9) +
  theme(strip.text.y.left = element_text(angle = 0, size = 7),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("FIG40_test1_rate_corrected.png", p1, width = 11, height = 9, dpi = 180, device = agg_png)

p2 <- NEPD %>% mutate(species = factor(species, SPORD)) %>%
  ggplot(aes(dS_nep, species, fill = species)) +
  geom_violin(scale = "width", colour = NA, alpha = .75) +
  geom_boxplot(width = .12, fill = "white", outlier.size = .3, linewidth = .3) +
  scale_y_discrete(limits = rev(SPORD)) +
  coord_cartesian(xlim = c(0, 3)) +
  labs(title = "The rate problem, shown directly",
       subtitle = "dS to Nepenthes. The true divergence time is IDENTICAL for all six, so any spread here is rate variation.",
       x = "dS to Nepenthes", y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "none")
ggsave("FIG41_rate_variation.png", p2, width = 9, height = 5, dpi = 180, device = agg_png)
cat("\nWROTE: FIG40_test1_rate_corrected.png FIG41_rate_variation.png\n")
