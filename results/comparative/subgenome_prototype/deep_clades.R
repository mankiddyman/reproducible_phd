#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Q1/Q2: does the deepest split within Droseraceae separate SUBGENOMES (each clade holding
# many species, each species spanning the split) or SPECIES (clades species-specific)?
# Also re-runs the per-species monophyly test on this exact tree set with a matched
# permutation null, to settle the binata full-tree-vs-quartet conflict.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ape); library(ggplot2); library(ragg); library(patchwork)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
set.seed(1); B <- 499
SPORD <- c("Dionaea_muscipula","Drosera_regia","Drosera_binata",
           "Drosera_paradoxa","Drosera_scorpioides","Drosera_capensis")

meta <- read_tsv("wgd7/tip_meta.tsv", show_col_types = FALSE)
key <- function(x) gsub("@", "_", gsub("['\"]", "", x))

prep <- function(f) {
  tr <- tryCatch(read.tree(f), error = function(e) NULL)
  if (is.null(tr) || Ntip(tr) < 6) return(NULL)
  m <- meta[match(key(tr$tip.label), key(meta$tip)), ]
  if (any(is.na(m$genome))) return(NULL)
  nep <- tr$tip.label[m$genome == "Nepenthes_gracilis"]
  if (length(nep) != 1) return(NULL)
  tr <- tryCatch(root(tr, nep, resolve.root = TRUE), error = function(e) NULL)
  if (is.null(tr)) return(NULL)
  tr <- drop.tip(tr, nep)
  m <- meta[match(key(tr$tip.label), key(meta$tip)), ]
  list(tr = tr, sp = m$genome, anchor = sub("\\.tre$", "", basename(f)))
}

## deepest split -> the two root clades
split2 <- function(tr) {
  k <- tr$edge[tr$edge[, 1] == Ntip(tr) + 1, 2]
  if (length(k) != 2) return(NULL)
  lapply(k, function(z) if (z <= Ntip(tr)) tr$tip.label[z] else extract.clade(tr, z)$tip.label)
}

## per-tree statistics, optionally on permuted species labels
stats_one <- function(tr, sp) {
  h <- split2(tr); if (is.null(h)) return(NULL)
  n1 <- length(h[[1]]); n2 <- length(h[[2]])
  if (min(n1, n2) < 2) return(NULL)                # a lone tip is not a deep split
  side <- ifelse(tr$tip.label %in% h[[1]], 1L, 2L)
  d <- tibble(species = sp, side = side)
  per <- d %>% group_by(species) %>%
    summarise(k = n(), in1 = sum(side == 1), in2 = sum(side == 2), .groups = "drop") %>%
    mutate(spans = in1 > 0 & in2 > 0, mono = (in1 == 0 | in2 == 0))
  tibble(n_sp_clade1 = n_distinct(d$species[d$side == 1]),
         n_sp_clade2 = n_distinct(d$species[d$side == 2]),
         n_sp = n_distinct(d$species),
         balance = min(n1, n2) / (n1 + n2),
         per = list(per))
}

fs <- list.files("wgd7/tre", "\\.tre$", full.names = TRUE)
P <- Filter(Negate(is.null), lapply(fs, prep))
cat(sprintf("trees usable: %d of %d\n", length(P), length(fs)))

OBS <- bind_rows(lapply(P, function(x) {
  s <- stats_one(x$tr, x$sp); if (is.null(s)) return(NULL)
  mutate(s, anchor = x$anchor) }))
cat(sprintf("trees with a balanced deep split: %d\n", nrow(OBS)))

## ---- Q1/Q2 headline: are the deep clades species-mixed? --------------------
cat("\n=== how many species does each deep clade contain? ===\n")
cat(sprintf("median species in the SMALLER clade: %.1f  (of %.1f present)\n",
            median(pmin(OBS$n_sp_clade1, OBS$n_sp_clade2)), median(OBS$n_sp)))
print(as.data.frame(OBS %>%
  mutate(minsp = pmin(n_sp_clade1, n_sp_clade2)) %>% count(minsp)))
cat("minsp = 1 means the deep split isolates ONE species -> species-specific event\n")
cat("minsp >= 3 means both clades hold many species -> shared subgenome structure\n")

## ---- per species: do its copies SPAN the deep split? -----------------------
PER <- bind_rows(lapply(seq_len(nrow(OBS)), function(i)
  mutate(OBS$per[[i]], anchor = OBS$anchor[i]))) %>% filter(k >= 2)

nullrun <- function() {
  bind_rows(lapply(P, function(x) {
    s <- stats_one(x$tr, sample(x$sp)); if (is.null(s)) return(NULL)
    s$per[[1]] })) %>% filter(k >= 2) %>%
    group_by(species) %>% summarise(f = mean(spans), .groups = "drop")
}
NUL <- bind_rows(lapply(seq_len(B), function(i) nullrun()))

R <- PER %>% group_by(species) %>%
  summarise(loci = n(), frac_spans = mean(spans), .groups = "drop") %>%
  left_join(NUL %>% group_by(species) %>%
              summarise(null_mean = mean(f), null_sd = sd(f), .groups = "drop"),
            by = "species") %>%
  mutate(z = (frac_spans - null_mean)/null_sd,
         verdict = ifelse(z > 2, "SPANS (shared A/B)",
                   ifelse(z < -2, "CLUSTERS (own event)", "at chance"))) %>%
  arrange(match(species, SPORD))
cat("\n=== do a species' own copies span the deepest Droseraceae split? ===\n")
print(as.data.frame(R %>% transmute(species, loci, frac_spans = round(frac_spans,3),
  null = round(null_mean,3), z = round(z,2), verdict)), row.names = FALSE)
cat("\nThis is the SAME question as the binata monophyly/quartet conflict,\n")
cat("run on one tree set with one null. Its answer supersedes both.\n")

## ---- pairwise: do two species co-occur in the same deep clade? -------------
CO <- bind_rows(lapply(seq_len(nrow(OBS)), function(i) {
  p <- OBS$per[[i]] %>% filter(k >= 1)
  if (nrow(p) < 2) return(NULL)
  cb <- combn(nrow(p), 2)
  tibble(a = p$species[cb[1,]], b = p$species[cb[2,]],
         both1 = p$in1[cb[1,]] > 0 & p$in1[cb[2,]] > 0,
         both2 = p$in2[cb[1,]] > 0 & p$in2[cb[2,]] > 0) })) %>%
  mutate(share = both1 | both2)
cat("\n=== fraction of trees where two species share at least one deep clade ===\n")
print(as.data.frame(CO %>% group_by(a, b) %>%
  summarise(n = n(), share = round(mean(share),3), .groups="drop") %>%
  arrange(desc(share))), row.names = FALSE)

write_csv(R, "deep_clade_spanning.csv"); write_csv(OBS %>% select(-per), "deep_clade_trees.csv")

write_csv(NUL, "deep_clade_null_reps.csv")

NUL2 <- NUL %>% mutate(species = factor(species, rev(SPORD)))
R2   <- R   %>% mutate(species = factor(species, rev(SPORD)))

pA <- ggplot() +
  geom_violin(data = NUL2, aes(f, species), fill = "grey80", colour = NA,
              scale = "width", width = .8) +
  geom_boxplot(data = NUL2, aes(f, species), width = .10, outlier.shape = NA,
               fill = "white", linewidth = .3) +
  geom_point(data = R2, aes(frac_spans, species, colour = verdict), size = 4.5) +
  geom_text(data = R2, aes(frac_spans, species,
            label = sprintf("z=%.1f", z)), vjust = -1.4, size = 3) +
  scale_colour_manual(values = c(`SPANS (shared A/B)` = "#1b9e77",
                                 `CLUSTERS (own event)` = "#d95f02",
                                 `at chance` = "grey45")) +
  labs(title = "Do a species' own gene copies span the deepest Droseraceae split?",
       subtitle = paste0("grey violin = ", B, " permutations of the species labels on that same tree; point = observed.\n",
                         "SPANS: copies sit on opposite sides, which shared A/B progenitors produce.\n",
                         "CLUSTERS: copies sit together, which a lineage-specific duplication produces."),
       x = "fraction of loci where the species' copies span the deep split",
       y = NULL, colour = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "top")

MS <- OBS %>% mutate(minsp = pmin(n_sp_clade1, n_sp_clade2)) %>% count(minsp) %>%
  mutate(frac = n/sum(n))
pB <- ggplot(MS, aes(factor(minsp), frac)) +
  geom_col(aes(fill = minsp >= 3), width = .7) +
  geom_text(aes(label = n), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c(`FALSE` = "#d95f02", `TRUE` = "#1b9e77"), guide = "none") +
  labs(title = "How many species does the SMALLER of the two deep clades contain?",
       subtitle = "1-2 (orange): the deep split mostly separates species  |  3+ (green): both clades are species-mixed, the shared-subgenome signature",
       x = "species in the smaller deep clade (of 6 possible)", y = "fraction of trees") +
  theme_bw(base_size = 10)

ggsave("FIG42_deep_clade_spanning.png", pA / pB + plot_layout(heights = c(1.5, 1)),
       width = 11, height = 10, dpi = 180, device = agg_png)
cat("\nWROTE: FIG42_deep_clade_spanning.png deep_clade_spanning.csv\n")
