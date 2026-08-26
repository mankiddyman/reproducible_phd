#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

d <- read_csv("distances.csv", show_col_types = FALSE) %>%
  mutate(a = sub("_sg[0-9]+_s[0-9]+$", "", chr1),
         b = sub("_sg[0-9]+_s[0-9]+$", "", chr2),
         class = ifelse(chr1 == chr2, "same_chr",
                 ifelse(a == b, "within_pair", "cross_pair")))

l <- bind_rows(
  transmute(d, class, comp = "Dionaea vs Dionaea (homeologs)", ks = dS_dio1_dio2),
  transmute(d, class, comp = "Dionaea vs Nepenthes", ks = dS_nep_dio1),
  transmute(d, class, comp = "Dionaea vs Nepenthes", ks = dS_nep_dio2)) %>%
  filter(is.finite(ks), ks > 0, ks < 3)

mode_of <- function(x) { dd <- density(x, from = 0, to = 2, n = 2048); dd$x[which.max(dd$y)] }
cat("=== all triplets ===\n")
print(as.data.frame(l %>% group_by(comp) %>%
  summarise(n = n(), median = median(ks), mode = mode_of(ks),
            q25 = quantile(ks, .25), q75 = quantile(ks, .75))), digits = 4)

cat("\n=== homeolog Ks by triplet class ===\n")
print(as.data.frame(l %>% filter(grepl("homeologs", comp)) %>% group_by(class) %>%
  summarise(n = n(), median = median(ks), mode = mode_of(ks),
            frac_lt_0.1 = mean(ks < 0.1))), digits = 4)

p1 <- ggplot(l, aes(ks, fill = comp, colour = comp)) +
  geom_histogram(binwidth = 0.02, position = "identity", alpha = .35, colour = NA) +
  geom_density(aes(y = after_stat(count) * 0.02), linewidth = .8, fill = NA) +
  coord_cartesian(xlim = c(0, 2)) +
  scale_fill_manual(values = c("#7570b3", "#d95f02")) +
  scale_colour_manual(values = c("#7570b3", "#d95f02")) +
  labs(title = "Ks distributions, Nepenthes-anchored triplets",
       subtitle = "binwidth 0.02; curves are scaled densities",
       x = "Ks (dS, NG86)", y = "gene pairs", fill = NULL, colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")

p2 <- l %>% filter(grepl("homeologs", comp)) %>%
  ggplot(aes(ks, fill = class)) +
  geom_histogram(binwidth = 0.02, position = "identity", alpha = .55, colour = NA) +
  facet_wrap(~ class, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(0, 1.5)) +
  scale_fill_manual(values = c(within_pair = "#7570b3",
                               cross_pair = "#d95f02", same_chr = "#1b9e77")) +
  labs(title = "Homeolog Ks split by triplet class",
       x = "Ks (dio1 vs dio2)", y = "gene pairs") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

ggsave("FIG11_ks_distributions.png", p1, width = 9, height = 5.5, dpi = 180, device = agg_png)
ggsave("FIG12_ks_by_class.png", p2, width = 8, height = 8, dpi = 180, device = agg_png)
cat("\nWROTE: FIG11_ks_distributions.png FIG12_ks_by_class.png\n")
