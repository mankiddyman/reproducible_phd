#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Is a species' own duplication OLDER or YOUNGER than its divergence from other species?
#   older  -> the duplication predates the speciations -> a shared ancient event
#   younger-> lineage-specific duplication
# Uses dS only. No trees, so the circularity in the span-vs-dS test does not apply.
# Rates largely cancel: dS(own)=2 r_X T_dup, dS(X,Y)=(r_X+r_Y) T_spec.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
SPORD <- c("Dionaea_muscipula","Drosera_regia","Drosera_binata",
           "Drosera_paradoxa","Drosera_scorpioides","Drosera_capensis")
NEP <- "Nepenthes_gracilis"

k <- read_csv("ks/pairwise_ks.csv", show_col_types = FALSE) %>%
  filter(!is.na(dS), dS > 0, dS < 5, codons >= 100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types = FALSE)
cn <- tm %>% count(nep_gene, genome, name = "copies")
md <- function(x, to = 3) { if (length(x) < 25) return(NA_real_)
  d <- density(x, from = 0, to = to, n = 2048); d$x[which.max(d$y)] }

## ---- speciation: SINGLE-COPY orthologs only, so no paralogy contaminates ----
sc <- cn %>% filter(copies == 1)
SPEC <- k %>% filter(sp1 != sp2, sp1 != NEP, sp2 != NEP) %>%
  semi_join(sc, by = c("anchor" = "nep_gene", "sp1" = "genome")) %>%
  semi_join(sc, by = c("anchor" = "nep_gene", "sp2" = "genome")) %>%
  mutate(a = pmin(sp1, sp2), b = pmax(sp1, sp2)) %>%
  group_by(a, b) %>%
  summarise(n = n(), spec_med = median(dS), spec_mode = md(dS), .groups = "drop")
cat("=== SPECIATION dS, single-copy orthologs only ===\n")
print(as.data.frame(SPEC %>% transmute(a, b, n, median = round(spec_med,3),
                                       mode = round(spec_mode,3))), row.names = FALSE)

## ---- duplication: within-species pairs -------------------------------------
DUP <- k %>% filter(sp1 == sp2, sp1 != NEP) %>% rename(species = sp1) %>%
  group_by(species) %>%
  summarise(n = n(), dup_med = median(dS), dup_mode = md(dS),
            dup_q25 = quantile(dS,.25), dup_q75 = quantile(dS,.75), .groups = "drop")
cat("\n=== DUPLICATION dS, a species' own copies ===\n")
print(as.data.frame(DUP %>% transmute(species, n, median = round(dup_med,3),
  mode = round(dup_mode,3), IQR = paste0(round(dup_q25,2),"-",round(dup_q75,2))) %>%
  arrange(match(species, SPORD))), row.names = FALSE)

## ---- the comparison ---------------------------------------------------------
CMP <- bind_rows(
  SPEC %>% transmute(species = a, other = b, spec_med),
  SPEC %>% transmute(species = b, other = a, spec_med)) %>%
  left_join(DUP %>% select(species, dup_med), by = "species") %>%
  mutate(ratio = dup_med / spec_med,
         verdict = ifelse(ratio > 1.15, "dup OLDER than split",
                   ifelse(ratio < 0.87, "dup YOUNGER than split", "~same")))
cat("\n=== duplication vs speciation, per species pair ===\n")
cat("ratio = median dS(own copies) / median dS(to that species)\n")
cat("  >1  duplication predates the split (shared / ancient)\n")
cat("  <1  duplication postdates the split (lineage-specific)\n\n")
print(as.data.frame(CMP %>% transmute(species, other,
  dup = round(dup_med,3), spec = round(spec_med,3),
  ratio = round(ratio,2), verdict) %>%
  arrange(match(species, SPORD), other)), row.names = FALSE)

cat("\n=== summary per species ===\n")
print(as.data.frame(CMP %>% group_by(species) %>%
  summarise(comparisons = n(), older = sum(ratio > 1.15),
            younger = sum(ratio < 0.87), median_ratio = round(median(ratio),2),
            .groups = "drop") %>% arrange(match(species, SPORD))), row.names = FALSE)
cat("\nA shared ancient A/B split predicts 'older' against EVERY other species.\n")
cat("Mixed results mean the duplication sits between some splits and not others,\n")
cat("which locates it on the phylogeny.\n")

## ---- per-locus version, so the spread is visible ---------------------------
# Per locus: the focal species has >=2 copies (that is the duplication) and the OTHER
# species is single-copy there (that is a clean speciation distance). Requiring the focal
# species to be single-copy too is contradictory, which is why the previous version was empty.
DUPL <- k %>% filter(sp1 == sp2, sp1 != NEP) %>% rename(species = sp1) %>%
  group_by(anchor, species) %>% summarise(dup = median(dS), .groups = "drop")

SPECL <- bind_rows(
    k %>% filter(sp1 != sp2, sp1 != NEP, sp2 != NEP) %>%
      transmute(anchor, species = sp1, other = sp2, spec = dS),
    k %>% filter(sp1 != sp2, sp1 != NEP, sp2 != NEP) %>%
      transmute(anchor, species = sp2, other = sp1, spec = dS)) %>%
  semi_join(sc, by = c("anchor" = "nep_gene", "other" = "genome")) %>%
  group_by(anchor, species, other) %>% summarise(spec = median(spec), .groups = "drop")

LOC <- inner_join(DUPL, SPECL, by = c("anchor", "species")) %>% mutate(older = dup > spec)
cat(sprintf("\nper-locus comparisons: %d\n", nrow(LOC)))

cat("\n=== per-locus: fraction where the duplication is OLDER than the split ===\n")
LT <- LOC %>% group_by(species) %>%
  summarise(loci = n(), n_older = sum(older), frac_older = round(mean(older),3),
            .groups = "drop") %>% arrange(match(species, SPORD)) %>% as.data.frame()
LT$p <- vapply(seq_len(nrow(LT)), function(i) {
  x <- as.integer(LT$n_older[i]); n <- as.integer(LT$loci[i])
  if (is.na(n) || n < 1 || is.na(x)) return(NA_real_)
  signif(binom.test(x, n, 0.5)$p.value, 3) }, numeric(1))
print(LT, row.names = FALSE)
cat("\nrows in LOC per species (0 means no single-copy comparison partner):\n")
print(as.data.frame(count(LOC, species)), row.names = FALSE)

p1 <- LOC %>% mutate(species = factor(species, SPORD)) %>%
  ggplot(aes(spec, dup)) +
  geom_abline(slope = 1, intercept = 0, colour = "#C0392B", linetype = 2) +
  geom_point(alpha = .18, size = .7, colour = "#5B4EA8") +
  facet_wrap(~ species, ncol = 3) +
  coord_cartesian(xlim = c(0,2.5), ylim = c(0,2.5)) +
  labs(title = "Is a species' own duplication older or younger than its speciation events?",
       subtitle = "one point per (locus, comparison species). Above the red line = duplication predates the split, which a shared ancient A/B event requires.",
       x = "dS to the other species at that locus (speciation)",
       y = "dS between this species' own copies (duplication)") +
  theme_bw(base_size = 9)
ggsave("FIG44_dup_vs_speciation.png", p1, width = 11, height = 7, dpi = 180, device = agg_png)
write_csv(CMP, "dup_vs_speciation.csv")
cat("\nWROTE: FIG44_dup_vs_speciation.png dup_vs_speciation.csv\n")
