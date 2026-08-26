#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Does the cross-lineage concordance survive when low-support votes are dropped?
# Votes with NA support (deciding node is the root or its child, so unlabelled) are
# KEPT -- their FastTree/IQ-TREE agreement is 0.953, i.e. they are not noise.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
set.seed(1); B <- 999

v0 <- read_csv("tract_votes7iq.csv", show_col_types = FALSE)

null_agree <- function(nl, pA, B = 999)
  mean(replicate(B, mean(sapply(nl, function(k) { x <- rbinom(1,k,pA); max(x,k-x)/k }))))

run <- function(thr) {
  v <- filter(v0, is.na(support) | support >= thr)
  lin <- v %>% group_by(pair, anchor, posA, sp, lineage) %>%
    summarise(fa = mean(vote == "A"), .groups = "drop") %>%
    filter(fa != 0.5) %>% mutate(call = ifelse(fa > .5, "A", "B"))
  g <- lin %>% group_by(pair, anchor) %>%
    summarise(nlin = n(), fA = mean(call == "A"),
              agree = max(fA, 1-fA), .groups = "drop") %>% filter(nlin >= 2)
  g %>% group_by(pair) %>%
    summarise(n = n(), mean_agree = mean(agree), pA = mean(fA),
              nl = list(nlin), .groups = "drop") %>%
    rowwise() %>% mutate(null = null_agree(nl[[1]], pA, B),
                         excess = mean_agree - null) %>% ungroup() %>%
    select(-nl) %>% mutate(thr = thr, votes = nrow(v))
}

R <- bind_rows(lapply(c(0, 50, 70, 95), run))
cat("=== cross-lineage concordance excess, by support threshold ===\n")
print(as.data.frame(R %>% select(thr, pair, n, mean_agree, null, excess) %>%
  pivot_wider(names_from = thr, values_from = c(n, excess),
              id_cols = pair)), digits = 3)
cat("\n=== votes retained and mean excess across pairs ===\n")
print(as.data.frame(R %>% group_by(thr) %>%
  summarise(votes = first(votes), anchors = sum(n),
            mean_excess = round(mean(excess), 4),
            min_excess = round(min(excess), 4), .groups = "drop")))
cat("\nexcess holding or rising as thr goes up => the signal is in the well-supported votes\n")
cat("excess falling => the apparent concordance was driven by unresolved nodes\n")
write_csv(R, "support_threshold_scan.csv")
