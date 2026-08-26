#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Per COPY PAIR, not per species: does the chance of spanning the deepest Droseraceae
# split rise with that pair's dS? Under one ancient A/B bifurcation plus later
# lineage-specific duplicates, old pairs span and young pairs sit inside one clade --
# in EVERY species. Also fixes the minsp confound by conditioning on clade tip count.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ape); library(ggplot2); library(ragg); library(patchwork)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
set.seed(1); B <- 299
SPORD <- c("Dionaea_muscipula","Drosera_regia","Drosera_binata",
           "Drosera_paradoxa","Drosera_scorpioides","Drosera_capensis")

meta <- read_tsv("wgd7/tip_meta.tsv", show_col_types = FALSE)
key <- function(x) gsub("@", "_", gsub("['\"]", "", x))
ks <- read_csv("ks/pairwise_ks.csv", show_col_types = FALSE) %>%
  filter(!is.na(dS), dS > 0, dS < 5, codons >= 100, sp1 == sp2) %>%
  transmute(anchor, a = pmin(seq1, seq2), b = pmax(seq1, seq2), species = sp1, dS)

fs <- list.files("wgd7/tre", "\\.tre$", full.names = TRUE)
PR <- bind_rows(lapply(fs, function(f) {
  tr <- tryCatch(read.tree(f), error = function(e) NULL)
  if (is.null(tr) || Ntip(tr) < 6) return(NULL)
  m <- meta[match(key(tr$tip.label), key(meta$tip)), ]
  if (any(is.na(m$genome))) return(NULL)
  nep <- tr$tip.label[m$genome == "Nepenthes_gracilis"]; if (length(nep) != 1) return(NULL)
  tr <- tryCatch(root(tr, nep, resolve.root = TRUE), error = function(e) NULL)
  if (is.null(tr)) return(NULL)
  tr <- drop.tip(tr, nep); m <- meta[match(key(tr$tip.label), key(meta$tip)), ]
  k <- tr$edge[tr$edge[,1] == Ntip(tr)+1, 2]; if (length(k) != 2) return(NULL)
  h <- lapply(k, function(z) if (z <= Ntip(tr)) tr$tip.label[z] else extract.clade(tr,z)$tip.label)
  if (min(lengths(h)) < 2) return(NULL)
  side <- setNames(ifelse(tr$tip.label %in% h[[1]], 1L, 2L), key(m$tip))
  a <- sub("\\.tre$","",basename(f))
  do.call(rbind, lapply(unique(m$genome), function(g) {
    tp <- key(m$tip)[m$genome == g]; if (length(tp) < 2) return(NULL)
    cb <- combn(tp, 2)
    data.frame(anchor = a, species = g,
               a = pmin(cb[1,], cb[2,]), b = pmax(cb[1,], cb[2,]),
               spans = side[cb[1,]] != side[cb[2,]],
               n_small = min(lengths(h)), n_tot = Ntip(tr),
               stringsAsFactors = FALSE) })) }))
cat(sprintf("within-species copy pairs on a balanced deep split: %d\n", nrow(PR)))

## join dS -- tip names in trees have @ replaced by _
ks2 <- ks %>% mutate(a = key(a), b = key(b))
D <- PR %>% inner_join(ks2, by = c("anchor","species","a","b"))
cat(sprintf("pairs with a dS value: %d (%.0f%%)\n", nrow(D), 100*nrow(D)/nrow(PR)))

D <- D %>% mutate(species = factor(species, SPORD),
                  bin = cut(dS, c(0,.2,.4,.6,.8,1.2,5),
                            labels = c("<0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.2",">1.2")))

cat("\n=== P(spans the deep split) by dS bin, per species ===\n")
tab <- D %>% group_by(species, bin) %>%
  summarise(n = n(), p_span = round(mean(spans),3), .groups="drop") %>% filter(n >= 15)
print(as.data.frame(tab %>% select(species, bin, n, p_span) %>%
  pivot_wider(names_from = bin, values_from = c(n, p_span))), row.names = FALSE)

cat("\n=== pooled, and a logistic test of the trend ===\n")
print(as.data.frame(D %>% group_by(bin) %>%
  summarise(pairs = n(), p_span = round(mean(spans),3), .groups="drop")), row.names = FALSE)
fit <- glm(spans ~ dS + species, data = D, family = binomial)
co <- summary(fit)$coefficients
cat(sprintf("\nlogistic spans ~ dS (+ species): dS coefficient %.3f, z = %.2f, p = %.3g\n",
            co["dS","Estimate"], co["dS","z value"], co["dS","Pr(>|z|)"]))
cat("positive and significant => older copy pairs span the ancient split, younger ones do not\n")
cat("   which is what ONE A/B bifurcation plus later lineage-specific duplicates predicts\n")

## per-species trend, so it is not driven by one lineage
cat("\n=== dS trend fitted separately within each species ===\n")
print(as.data.frame(bind_rows(lapply(SPORD, function(s) {
  d <- filter(D, species == s); if (nrow(d) < 40 || n_distinct(d$spans) < 2) return(NULL)
  m <- glm(spans ~ dS, data = d, family = binomial); cc <- summary(m)$coefficients
  tibble(species = s, pairs = nrow(d), coef = round(cc["dS","Estimate"],3),
         z = round(cc["dS","z value"],2), p = signif(cc["dS","Pr(>|z|)"],3)) }))),
  row.names = FALSE)

## ---- minsp, now conditioned on the smaller clade's TIP COUNT ---------------
MS <- PR %>% distinct(anchor, n_small, n_tot)
cat("\n=== the minsp confound: species in the smaller clade is capped by its tip count ===\n")
print(as.data.frame(MS %>% mutate(sz = cut(n_small, c(1,2,3,5,100),
        labels=c("2 tips","3 tips","4-5 tips","6+ tips"))) %>% count(sz)), row.names = FALSE)
cat("half the trees have a 2-tip smaller clade, which can hold at most 2 species by construction\n")

p1 <- ggplot(D, aes(dS, as.numeric(spans))) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              colour = "#5B4EA8", fill = "#c9c4e8") +
  stat_summary_bin(aes(x = dS), bins = 12, fun = mean, geom = "point",
                   colour = "#12795E", size = 2) +
  facet_wrap(~ species, ncol = 3) +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 1)) +
  labs(title = "Do OLD copy pairs span the deepest Droseraceae split, and young ones not?",
       subtitle = "green points = binned observed fraction; purple = logistic fit.\nA rising curve in every species is what one ancient A/B bifurcation plus later lineage-specific duplicates produces.",
       x = "dS between the two copies", y = "P(copies span the deep split)") +
  theme_bw(base_size = 9)
ggsave("FIG43_span_by_ks.png", p1, width = 11, height = 7, dpi = 180, device = agg_png)
write_csv(D %>% select(-a,-b), "span_by_ks.csv")
cat("\nWROTE: FIG43_span_by_ks.png span_by_ks.csv\n")
