#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# z = (obs - null_mu)/null_sd, computed without rowwise list-columns.
# Downsample control: same number of votes, drawn at random from the unfiltered set.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
set.seed(1); B <- 999; NDOWN <- 20

v0 <- read_csv("tract_votes7iq.csv", show_col_types = FALSE)

stat <- function(v) {
  lin <- v %>% group_by(pair, anchor, sp, lineage) %>%
    summarise(fa = mean(vote == "A"), .groups = "drop") %>%
    filter(fa != 0.5) %>% mutate(call = ifelse(fa > .5, 1L, 0L))
  g <- lin %>% group_by(pair, anchor) %>%
    summarise(nlin = n(), x = sum(call), .groups = "drop") %>% filter(nlin >= 2) %>%
    mutate(agree = pmax(x, nlin - x)/nlin)
  if (!nrow(g)) return(NULL)
  do.call(rbind, lapply(split(g, g$pair), function(d) {
    pA <- sum(d$x)/sum(d$nlin)
    nulls <- vapply(seq_len(B), function(i) {
      xx <- rbinom(nrow(d), d$nlin, pA)
      mean(pmax(xx, d$nlin - xx)/d$nlin) }, 0)
    obs <- mean(d$agree); mu <- mean(nulls); sdv <- sd(nulls)
    data.frame(pair = d$pair[1], n = nrow(d), obs = obs, null_mu = mu, null_sd = sdv,
               excess = obs - mu, z = if (sdv > 0) (obs - mu)/sdv else NA_real_,
               p = (1 + sum(nulls >= obs))/(B + 1), stringsAsFactors = FALSE)
  }))
}

out <- bind_rows(lapply(c(0, 50, 70, 95), function(thr) {
  vf <- filter(v0, is.na(support) | support >= thr)
  s  <- stat(vf); if (is.null(s)) return(NULL)
  s$thr <- thr; s$kind <- "filtered"; s$votes <- nrow(vf)
  dn <- bind_rows(lapply(seq_len(NDOWN), function(i) stat(slice_sample(v0, n = nrow(vf)))))
  ds <- dn %>% group_by(pair) %>%
    summarise(n = mean(n), obs = mean(obs), null_mu = mean(null_mu),
              null_sd = mean(null_sd), excess = mean(excess),
              z = mean(z, na.rm = TRUE), p = mean(p), .groups = "drop") %>%
    mutate(thr = thr, kind = "downsample", votes = nrow(vf))
  bind_rows(s, ds)
}))
write_csv(out, "support_fair.csv")

cat("=== z per pair: filtered vs size-matched random ===\n")
print(as.data.frame(out %>% select(pair, thr, kind, z) %>%
  pivot_wider(names_from = c(kind, thr), values_from = z)), digits = 3)

cat("\n=== usable anchors: filtered vs random of the same vote count ===\n")
print(as.data.frame(out %>% select(pair, thr, kind, n) %>%
  pivot_wider(names_from = c(kind, thr), values_from = n)), digits = 3)

cat("\n=== summary ===\n")
print(as.data.frame(out %>% group_by(thr, kind) %>%
  summarise(votes = first(votes), anchors = round(sum(n)),
            mean_z = round(mean(z, na.rm = TRUE), 2),
            min_z = round(min(z, na.rm = TRUE), 2),
            mean_excess = round(mean(excess), 4),
            sig = sum(p < 0.05), .groups = "drop") %>%
  arrange(thr, desc(kind))), digits = 3)

d <- out %>% select(pair, thr, kind, z) %>%
  pivot_wider(names_from = kind, values_from = z) %>%
  mutate(gain = filtered - downsample)
cat("\n=== filtered minus downsample z (positive = the filter adds information) ===\n")
print(as.data.frame(d %>% group_by(thr) %>%
  summarise(mean_gain = round(mean(gain, na.rm = TRUE), 3),
            pairs_better = sum(gain > 0, na.rm = TRUE), n = n(), .groups = "drop")))
