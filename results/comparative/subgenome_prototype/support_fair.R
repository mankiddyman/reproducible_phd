#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Excess is not comparable across thresholds because filtering changes the number of
# lineages per anchor, which moves the null. Two fixes:
#  (a) z-score: excess divided by the SD of the permuted null -> comparable
#  (b) DOWNSAMPLE control: at each threshold, also draw a random subset of the SAME
#      size from the UNFILTERED votes. If the filter beats its own downsample, the
#      support filter is doing real work rather than just shrinking n.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
set.seed(1); B <- 999; NDOWN <- 20

v0 <- read_csv("tract_votes7iq.csv", show_col_types = FALSE)

stat <- function(v) {
  lin <- v %>% group_by(pair, anchor, sp, lineage) %>%
    summarise(fa = mean(vote == "A"), .groups = "drop") %>%
    filter(fa != 0.5) %>% mutate(call = ifelse(fa > .5, "A", "B"))
  g <- lin %>% group_by(pair, anchor) %>%
    summarise(nlin = n(), fA = mean(call == "A"),
              agree = max(fA, 1-fA), .groups = "drop") %>% filter(nlin >= 2)
  if (!nrow(g)) return(NULL)
  g %>% group_by(pair) %>%
    summarise(n = n(), obs = mean(agree), pA = mean(fA), nl = list(nlin), .groups="drop") %>%
    rowwise() %>%
    mutate(nulls = list(replicate(B, mean(sapply(nl[[1]], function(k) {
             x <- rbinom(1,k,pA); max(x,k-x)/k })))),
           null_mu = mean(nulls[[1]]), null_sd = sd(nulls[[1]]),
           excess = obs - null_mu,
           z = (obs - null_mu)/null_sd,
           p = (1 + sum(nulls[[1]] >= obs))/(B+1)) %>%
    ungroup() %>% select(pair, n, obs, null_mu, null_sd, excess, z, p)
}

out <- bind_rows(lapply(c(0, 50, 70, 95), function(thr) {
  vf <- filter(v0, is.na(support) | support >= thr)
  s  <- stat(vf); if (is.null(s)) return(NULL)
  s$thr <- thr; s$kind <- "filtered"; s$votes <- nrow(vf)
  # size-matched random control from the unfiltered set
  dn <- bind_rows(lapply(seq_len(NDOWN), function(i) {
    sd_ <- stat(slice_sample(v0, n = nrow(vf)))
    if (is.null(sd_)) NULL else mutate(sd_, rep = i) }))
  ds <- dn %>% group_by(pair) %>%
    summarise(n = round(mean(n)), obs = mean(obs), null_mu = mean(null_mu),
              null_sd = mean(null_sd), excess = mean(excess), z = mean(z),
              p = mean(p), .groups="drop") %>%
    mutate(thr = thr, kind = "downsample", votes = nrow(vf))
  bind_rows(s, ds)
}))
write_csv(out, "support_fair.csv")

cat("=== z-score (excess / null SD) — comparable across thresholds ===\n")
print(as.data.frame(out %>% select(pair, thr, kind, z) %>%
  pivot_wider(names_from = c(kind, thr), values_from = z)), digits = 3)

cat("\n=== summary: filter vs a size-matched random subset ===\n")
print(as.data.frame(out %>% group_by(thr, kind) %>%
  summarise(votes = first(votes), anchors = sum(n),
            mean_z = round(mean(z),2), min_z = round(min(z),2),
            pairs_p_lt_.05 = sum(p < 0.05), .groups="drop") %>%
  arrange(thr, desc(kind))), digits = 3)

cat("\nfiltered z > downsample z at the same thr => the support filter adds information\n")
cat("filtered z ~ downsample z              => it is only reducing n\n")
