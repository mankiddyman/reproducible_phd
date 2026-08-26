#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Per-species subgenome-affinity profiles, the distances between them, and a
# bootstrapped clustering. Agreement between two species is invariant to the arbitrary
# sgA/sgB flip within a Dionaea pair, so distances are well defined; the ORIENTATION
# (which side is which progenitor) is not, and is never used here.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ape); library(patchwork)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
set.seed(1); NBOOT <- 1000; MINSHARE <- 30
SPORD <- c("regia","binata","paradoxa","scorpioides","capensis")

v <- read_csv("tract_votes_blocks7.csv", show_col_types = FALSE) %>%
  mutate(species = sub("Drosera_","",sp))

## majority call per (species, anchor); ties dropped
call <- v %>% group_by(pair, anchor, species) %>%
  summarise(k = n(), a = sum(vote == "A"), .groups = "drop") %>%
  filter(a != k - a) %>% mutate(call = ifelse(a > k - a, 1L, 0L))
cat(sprintf("votes %d | (species, anchor) calls %d | anchors %d\n",
            nrow(v), nrow(call), n_distinct(call$anchor)))

## ---- 1. summary per species per Dionaea pair -------------------------------
S <- call %>% group_by(species, pair) %>%
  summarise(genes = n(), fracA = mean(call), .groups = "drop")
Sall <- call %>% group_by(species) %>%
  summarise(genes = n(), fracA = round(mean(call), 3),
            copies_per_gene = round(nrow(filter(v, species == species[1])) / n(), 2),
            .groups = "drop") %>% arrange(match(species, SPORD))
cat("\n=== overall affinity per species (fracA is direction-arbitrary per pair) ===\n")
print(as.data.frame(Sall))
cat("\n=== fraction voting sgA, species x Dionaea pair ===\n")
print(as.data.frame(S %>% mutate(fracA = round(fracA,2)) %>%
  select(species, pair, fracA) %>% pivot_wider(names_from = pair, values_from = fracA) %>%
  arrange(match(species, SPORD))))
write_csv(S, "species_pair_profiles.csv")

## ---- 2. pairwise agreement, and excess over independence -------------------
W <- call %>% select(pair, anchor, species, call) %>%
  pivot_wider(names_from = species, values_from = call)
sp <- intersect(SPORD, names(W))
agr <- exc <- shr <- matrix(NA_real_, length(sp), length(sp), dimnames = list(sp, sp))
for (i in seq_along(sp)) for (j in seq_along(sp)) {
  x <- W[[sp[i]]]; y <- W[[sp[j]]]; ok <- !is.na(x) & !is.na(y)
  shr[i,j] <- sum(ok)
  if (sum(ok) >= MINSHARE) {
    a <- mean(x[ok] == y[ok]); px <- mean(x[ok]); py <- mean(y[ok])
    agr[i,j] <- a; exc[i,j] <- a - (px*py + (1-px)*(1-py))
  }
}
cat("\n=== shared anchor genes between species ===\n"); print(shr)
cat("\n=== observed agreement ===\n"); print(round(agr, 3))
cat("\n=== EXCESS agreement over independence (the informative quantity) ===\n")
print(round(exc, 3))
write.csv(round(exc,4), "species_excess_agreement.csv")

## ---- 3. clustering with bootstrap over anchor genes ------------------------
dmat <- function(Wi) {
  m <- matrix(0, length(sp), length(sp), dimnames = list(sp, sp))
  for (i in seq_along(sp)) for (j in seq_along(sp)) {
    x <- Wi[[sp[i]]]; y <- Wi[[sp[j]]]; ok <- !is.na(x) & !is.na(y)
    m[i,j] <- if (sum(ok) >= 10) 1 - mean(x[ok] == y[ok]) else NA_real_ }
  m[is.na(m)] <- mean(m, na.rm = TRUE); m
}
D <- dmat(W); hc <- hclust(as.dist(D), method = "average"); tr <- as.phylo(hc)
bips <- function(t) sapply(prop.part(t), function(z) paste(sort(t$tip.label[z]), collapse=","))
obs <- bips(tr)
hits <- setNames(rep(0L, length(obs)), obs)
for (b in seq_len(NBOOT)) {
  Wb <- W[sample(nrow(W), nrow(W), replace = TRUE), ]
  tb <- tryCatch(as.phylo(hclust(as.dist(dmat(Wb)), method="average")), error=function(e) NULL)
  if (is.null(tb)) next
  h <- bips(tb); for (k in names(hits)) if (k %in% h) hits[k] <- hits[k] + 1L
}
cat(sprintf("\n=== clustering by subgenome-affinity profile (%d bootstraps) ===\n", NBOOT))
print(data.frame(clade = names(hits), support = round(100*hits/NBOOT)), row.names = FALSE)
cat("\nNewick: "); cat(write.tree(tr)); cat("\n")
cat("NOTE: this clusters by which homeolog each species retained, not by phylogeny.\n")

## ---- 4. figure --------------------------------------------------------------
p1 <- S %>% mutate(species = factor(species, rev(SPORD))) %>%
  ggplot(aes(pair, species, fill = fracA)) +
  geom_tile(colour = "white", linewidth = .5) +
  geom_text(aes(label = sprintf("%.2f\n(%d)", fracA, genes)), size = 2.6) +
  scale_fill_gradient2(low = "#d95f02", mid = "grey93", high = "#1b9e77",
                       midpoint = .5, limits = c(0,1)) +
  labs(title = "Fraction of genes placed with the local sgA",
       subtitle = "sgA/sgB are arbitrary WITHIN a column; compare down columns, not across",
       x = "Dionaea chromosome pair", y = NULL, fill = "frac sgA") +
  theme_bw(base_size = 9)

ex <- as.data.frame(as.table(exc)); names(ex) <- c("a","b","excess")
p2 <- ex %>% filter(!is.na(excess)) %>%
  mutate(a = factor(a, SPORD), b = factor(b, rev(SPORD))) %>%
  ggplot(aes(a, b, fill = excess)) +
  geom_tile(colour="white", linewidth=.5) +
  geom_text(aes(label = sprintf("%.3f", excess)), size = 2.8) +
  scale_fill_gradient2(low="#5B4EA8", mid="grey95", high="#B7791F", midpoint=0) +
  labs(title = "Excess agreement over independence",
       subtitle = "how much more two species retain the same homeolog than their marginal rates predict",
       x = NULL, y = NULL, fill = "excess") +
  theme_bw(base_size = 9)

png("FIG32_species_profiles.png", width = 2400, height = 1500, res = 190)
print(p1 / p2); dev.off()

pdf("FIG33_species_dendrogram.pdf", width = 8, height = 5.5)
plot(tr, cex = 1.1, edge.width = 2, main = "Clustering by subgenome-affinity profile")
nodelabels(sprintf("%d", round(100*hits[obs])), frame = "none", adj = c(1.2, -0.4), cex = .85)
add.scale.bar()
mtext("distance = fraction of shared anchor genes with a different homeolog call; numbers = bootstrap %",
      side = 1, line = 3, cex = .75)
dev.off()
cat("\nWROTE: FIG32_species_profiles.png FIG33_species_dendrogram.pdf\n")
