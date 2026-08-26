#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(patchwork); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
th <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())

rc   <- read_csv("raw_gene_content.csv", show_col_types = FALSE)
bl   <- read_csv("blast_presence_summary.csv", show_col_types = FALSE)
v    <- read_csv("topology_votes.csv", show_col_types = FALSE)
pw   <- read_csv("topology_rate_corrected.csv", show_col_types = FALSE)
frac <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE)
ds   <- read_csv("phasing_by_chrpair_v2.csv", show_col_types = FALSE) %>% filter(set == "conv_filtered")

part <- frac %>% filter(anchor == "Nepenthes_gracilis") %>%
  transmute(pair_lab, X = winner, Y = ifelse(winner == chrA, chrB, chrA))

## A. fractionation ratio under stricter filters
pA <- rc %>% select(pair_lab, r_all, r_keep, r_cons, r_hog) %>%
  pivot_longer(-pair_lab, names_to = "f", values_to = "ratio") %>%
  mutate(f = factor(f, c("r_all","r_keep","r_cons","r_hog"),
                    c("all genes","OMArk keep","OMArk consistent","in orthogroup"))) %>%
  ggplot(aes(f, ratio, group = pair_lab, colour = pair_lab)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  geom_line(alpha = .8) + geom_point(size = 1.6) +
  scale_colour_brewer(palette = "Dark2") +
  labs(title = "A. Gene-count ratio X/Y survives every quality filter",
       subtitle = "Stricter filters raise the ratio — an annotation artifact would lower it",
       x = NULL, y = "genes on X / genes on Y", colour = NULL) +
  th + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "right")

## B. tblastn presence
pB <- bl %>%
  mutate(cat = factor(cat, c("posX","posY","lostX","lostY","negX"),
                      c("control:\npresent on X","control:\npresent on Y",
                        "X lost it","Y lost it","background")),
         grp = ifelse(grepl("control", cat), "positive control", "test")) %>%
  ggplot(aes(cat, frac_present, fill = grp)) +
  geom_boxplot(outlier.shape = NA, alpha = .5, width = .6) +
  geom_jitter(width = .12, size = 1.4, alpha = .8) +
  scale_fill_manual(values = c("positive control" = "#1b9e77", "test" = "#d95f02")) +
  ylim(0, 1) +
  labs(title = "B. Missing genes are physically absent, not unannotated",
       subtitle = "each point = one chromosome pair",
       x = NULL, y = "fraction with intact match", fill = NULL) +
  th + theme(legend.position = "none")

## C. topology marginal skew, Nepenthes as control
pC <- v %>% group_by(sp) %>%
  summarise(n = n(), k = sum(rel > 0), .groups = "drop") %>%
  rowwise() %>% mutate(f = k/n, lo = binom.test(k, n)$conf.int[1],
                       hi = binom.test(k, n)$conf.int[2]) %>% ungroup() %>%
  mutate(ctrl = ifelse(sp == "Nepenthes_gracilis", "control (outside both)", "Drosera")) %>%
  ggplot(aes(f, reorder(sp, f), colour = ctrl)) +
  geom_vline(xintercept = .5, linetype = 2, colour = "grey40") +
  geom_linerange(aes(xmin = lo, xmax = hi)) + geom_point(aes(size = n)) +
  scale_colour_manual(values = c("Drosera" = "#7570b3", "control (outside both)" = "#666666")) +
  scale_size_continuous(range = c(1.5, 3.5), guide = "none") +
  labs(title = "C. Every Drosera prefers X; Nepenthes does not",
       subtitle = "fraction of genes scoring higher against X (95% CI)",
       x = "fraction voting X", y = NULL, colour = NULL) +
  th + theme(legend.position = "bottom")

## D. rate-corrected topology per pair
hl <- function(x) { w <- suppressWarnings(wilcox.test(x, conf.int = TRUE))
  tibble(est = unname(w$estimate), lo = w$conf.int[1], hi = w$conf.int[2]) }
pD <- pw %>% group_by(pair_lab) %>% reframe(n = n(), hl(delta)) %>%
  ggplot(aes(est, reorder(pair_lab, est))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_linerange(aes(xmin = lo, xmax = hi), colour = "#7570b3") +
  geom_point(aes(size = n), colour = "#7570b3") +
  scale_size_continuous(range = c(1.5, 3.5)) +
  labs(title = "D. Rate-corrected topology signal (Drosera minus Nepenthes)",
       subtitle = "positive = Drosera closer to X than rate alone explains",
       x = "median per-HOG delta", y = NULL, size = "HOGs") + th

## E. sign vector across methods
sub_blue <- c("chr1_sg1_s5","chr8_sg2_s14","chr6_sg1_s6",
              "chr4_sg2_s12","chr5_sg2_s16","chr7_sg2_s13")
topo <- pw %>% group_by(pair_lab) %>% summarise(d = median(delta), .groups = "drop")
sv <- part %>%
  left_join(select(ds, pair_lab, faster), by = "pair_lab") %>%
  left_join(topo, by = "pair_lab") %>%
  mutate(fractionation = "X",
         `dS (faster)`  = ifelse(faster == X, "X", "Y"),
         topology       = ifelse(d > 0, "X", "Y"),
         `SubPhaser`    = case_when(X %in% sub_blue ~ "X", Y %in% sub_blue ~ "Y",
                                    TRUE ~ NA_character_)) %>%
  select(pair_lab, fractionation, `dS (faster)`, topology, SubPhaser) %>%
  pivot_longer(-pair_lab, names_to = "method", values_to = "pick") %>%
  mutate(method = factor(method, c("fractionation","topology","dS (faster)","SubPhaser")))
pE <- ggplot(sv, aes(method, pair_lab, fill = pick)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(is.na(pick), "unresolved", pick)), size = 3) +
  scale_fill_manual(values = c(X = "#a6dba0", Y = "#f4a582"), na.value = "grey88") +
  labs(title = "E. Which member does each method pick?",
       subtitle = "X is defined by fractionation, so column 1 is X by construction",
       x = NULL, y = NULL, fill = NULL) + th + theme(legend.position = "none")

fig1 <- (pA | pB) / (pC | pD) / pE + plot_layout(heights = c(1, 1, .9))
ggsave("FIG1_evidence_summary.png", fig1, width = 13, height = 13, dpi = 150, device = agg_png)
ggsave("FIG1_evidence_summary.pdf", fig1, width = 13, height = 13)

## ---- FIG 2: topology detail, H2 vs H3 --------------------------------------
pF <- v %>% mutate(grp = ifelse(sp == "Nepenthes_gracilis", "Nepenthes (control)", "Drosera")) %>%
  ggplot(aes(rel, colour = sp, linetype = grp)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_density(linewidth = .7) + coord_cartesian(xlim = c(-.12, .12)) +
  scale_colour_brewer(palette = "Dark2") +
  labs(title = "F. Distribution of per-gene preference",
       subtitle = "unimodal shift = one lineage closer to X; bimodal = shared WGD",
       x = "relative bitscore (bX - bY) / max", y = "density",
       colour = NULL, linetype = NULL) + th

## within-HOG concordance among multi-copy Drosera genes
mc <- v %>% filter(grepl("^Drosera", sp)) %>%
  group_by(sp, globHOG, pair_lab) %>% filter(n() >= 2) %>%
  summarise(n = n(), fx = mean(rel > 0), .groups = "drop") %>%
  mutate(concordant = fx == 0 | fx == 1)
cc <- mc %>% group_by(sp) %>%
  summarise(n_hog = n(), frac_conc = mean(concordant), .groups = "drop")
pG <- ggplot(cc, aes(frac_conc, reorder(sp, frac_conc))) +
  geom_vline(xintercept = .5, linetype = 2, colour = "grey40") +
  geom_col(fill = "#7570b3", width = .6) +
  geom_text(aes(label = sprintf("n=%d", n_hog)), hjust = -0.15, size = 3) +
  xlim(0, 1) +
  labs(title = "G. Do multiple copies in one species agree?",
       subtitle = "high = nested lineage (H3); ~0.5 or lower = shared WGD (H2)",
       x = "fraction of HOGs where all copies vote the same way", y = NULL) + th

pH <- ggplot(pw, aes(delta)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_histogram(bins = 60, fill = "#1b9e77", colour = NA) +
  facet_wrap(~ pair_lab, ncol = 4, scales = "free_y") +
  coord_cartesian(xlim = c(-.1, .1)) +
  labs(title = "H. Per-HOG rate-corrected signal, by chromosome pair",
       x = "delta (Drosera - Nepenthes)", y = "HOGs") + th

fig2 <- (pF | pG) / pH + plot_layout(heights = c(1, 1))
ggsave("FIG2_topology_detail.png", fig2, width = 13, height = 9, dpi = 150, device = agg_png)
ggsave("FIG2_topology_detail.pdf", fig2, width = 13, height = 9)

cat("\n=== H2 vs H3: within-HOG concordance of multi-copy Drosera ===\n")
print(as.data.frame(cc), digits = 3)
cat("\nexpected under random voting: ~0.5 for 2-copy HOGs\n")
two <- mc %>% filter(n == 2) %>% group_by(sp) %>%
  summarise(n = n(), k = sum(concordant),
            p = binom.test(sum(concordant), n(), 0.5)$p.value, .groups = "drop") %>%
  mutate(frac = k/n)
cat("\n=== restricted to HOGs with exactly 2 copies ===\n")
print(as.data.frame(two), digits = 3)

cat("\nWROTE: FIG1_evidence_summary.{png,pdf}  FIG2_topology_detail.{png,pdf}\n")
