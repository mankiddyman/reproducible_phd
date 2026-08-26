#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# TEST 1: within-species dS. Split by how many copies that species has at the locus,
# because a k-copy locus contributes k(k-1)/2 pairwise comparisons and mixing them
# collapses several real modes into one meaningless peak.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ragg); library(patchwork)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
SPORD <- c("Dionaea_muscipula","Drosera_regia","Drosera_binata",
           "Drosera_paradoxa","Drosera_scorpioides","Drosera_capensis")
SAT <- 1.0   # above this, dS is saturated and unreliable

k <- read_csv("ks/pairwise_ks.csv", show_col_types = FALSE)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types = FALSE)
cn <- tm %>% count(nep_gene, genome, name = "copies")

cat(sprintf("raw pairs %d | NA dS %d (%.1f%%) | dS<=0 or >=5 %d (%.1f%%)\n",
            nrow(k), sum(is.na(k$dS)), 100*mean(is.na(k$dS)),
            sum(k$dS<=0 | k$dS>=5, na.rm=TRUE), 100*mean(k$dS<=0 | k$dS>=5, na.rm=TRUE)))

W <- k %>% filter(!is.na(dS), dS > 0, dS < 5, codons >= 100, sp1 == sp2) %>%
  rename(species = sp1) %>%
  left_join(cn, by = c("anchor" = "nep_gene", "species" = "genome")) %>%
  filter(!is.na(copies)) %>%
  mutate(kgrp = ifelse(copies >= 4, "4+ copies", paste0(copies, " copies")),
         species = factor(species, SPORD))
cat(sprintf("usable within-species pairs: %d\n", nrow(W)))

cat("\n=== pairs per species per copy number ===\n")
print(as.data.frame(W %>% count(species, kgrp) %>%
  pivot_wider(names_from = kgrp, values_from = n, values_fill = 0)))

md <- function(x, from=0, to=2) { if (length(x) < 20) return(NA_real_)
  d <- density(x, from=from, to=to, n=2048); d$x[which.max(d$y)] }
cat("\n=== dS summary, split by copy number ===\n")
print(as.data.frame(W %>% group_by(species, kgrp) %>%
  summarise(n = n(), mode = round(md(dS),3), median = round(median(dS),3),
            frac_saturated = round(mean(dS > SAT),2), .groups="drop") %>%
  filter(n >= 20)), row.names = FALSE)

dio <- md(filter(W, species=="Dionaea_muscipula", copies==2)$dS)
cat(sprintf("\nDionaea 2-copy reference mode: %.3f\n", dio))

p1 <- ggplot(W, aes(dS)) +
  annotate("rect", xmin=SAT, xmax=2, ymin=-Inf, ymax=Inf, fill="grey85", alpha=.6) +
  geom_histogram(aes(y = after_stat(density)), binwidth=.05, boundary=0,
                 fill="#7570b3", colour=NA, alpha=.55) +
  geom_density(colour="#5B4EA8", linewidth=.7, bw=.06) +
  geom_vline(xintercept = dio, colour="#C0392B", linetype=2, linewidth=.6) +
  facet_grid(species ~ kgrp, scales="free_y", switch="y") +
  coord_cartesian(xlim=c(0, 2)) +
  labs(title="TEST 1 — divergence between a species' own gene copies (dS)",
       subtitle=paste0("dashed red = the Dionaea 2-copy mode (", round(dio,3),
                       "), i.e. when the two Dionaea progenitors split.\n",
                       "Grey zone = dS > 1, where synonymous sites saturate and estimates are lower bounds.\n",
                       "Columns split by how many copies that species has at the locus, because a k-copy locus gives k(k-1)/2 comparisons."),
       x="dS (synonymous divergence)", y=NULL) +
  theme_bw(base_size=9) +
  theme(strip.text.y.left = element_text(angle=0, size=7),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("FIG38_test1_within_species_dS.png", p1, width=11, height=9, dpi=180, device=agg_png)

p2 <- W %>% filter(dS < 2) %>%
  ggplot(aes(dS, species, fill = species)) +
  annotate("rect", xmin=SAT, xmax=2, ymin=-Inf, ymax=Inf, fill="grey85", alpha=.6) +
  geom_violin(scale="width", colour=NA, alpha=.75) +
  geom_boxplot(width=.12, outlier.size=.3, fill="white", linewidth=.3) +
  geom_vline(xintercept=dio, colour="#C0392B", linetype=2) +
  scale_y_discrete(limits=rev(SPORD)) +
  labs(title="All within-species comparisons pooled",
       subtitle="a single peak per species is only meaningful if that species has 2 copies",
       x="dS", y=NULL) +
  theme_bw(base_size=9) + theme(legend.position="none")
ggsave("FIG39_test1_violin.png", p2, width=9, height=5, dpi=180, device=agg_png)
cat("\nWROTE: FIG38_test1_within_species_dS.png FIG39_test1_violin.png\n")
