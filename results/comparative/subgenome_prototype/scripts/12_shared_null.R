#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); NPERM <- 999
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"

k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS > 0, dS < 5, codons >= 100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
cn <- tm %>% count(nep_gene, genome, name="copies")
sc <- cn %>% filter(copies == 1)

hr("0. logic")
cat("  Shared A/B: X's A and Y's A are ORTHOLOGS, so the smallest cross-species\n")
cat("  dS at a locus is a speciation distance. Independent: the donor sublineages\n")
cat("  had already split, so it exceeds speciation.\n")
cat("  BUT min-of-m draws sits below the median of the parent distribution.\n")
cat("  The null below generates that bias explicitly instead of assuming ratio 1.\n")

## empirical speciation distribution per species pair (single-copy orthologs)
SP <- k %>% filter(sp1 != sp2, sp1 != NEP, sp2 != NEP) %>%
  semi_join(sc, by=c("anchor"="nep_gene","sp1"="genome")) %>%
  semi_join(sc, by=c("anchor"="nep_gene","sp2"="genome")) %>%
  mutate(a=pmin(sp1,sp2), b=pmax(sp1,sp2)) %>% select(a, b, dS)

## observed: per locus, min cross-species dS where BOTH have >=2 copies
OB <- k %>% filter(sp1 != sp2, sp1 != NEP, sp2 != NEP) %>%
  mutate(a=pmin(sp1,sp2), b=pmax(sp1,sp2)) %>%
  inner_join(cn %>% rename(a=genome, ka=copies), by=c("anchor"="nep_gene","a")) %>%
  inner_join(cn %>% rename(b=genome, kb=copies), by=c("anchor"="nep_gene","b")) %>%
  filter(ka >= 2, kb >= 2) %>%
  group_by(a, b, anchor, ka, kb) %>%
  summarise(obs_min = min(dS), .groups="drop") %>%
  mutate(m = pmin(ka, kb))

hr("1. observed vs simulated-null minimum")
cat("  Null draws m values from that pair's own speciation distribution and\n")
cat("  takes the minimum, m = min(copies). That IS the shared-A/B prediction.\n\n")
set.seed(1)
out <- bind_rows(lapply(split(OB, paste(OB$a, OB$b)), function(d) {
  pool <- SP$dS[SP$a == d$a[1] & SP$b == d$b[1]]
  if (length(pool) < 30 || nrow(d) < 20) return(NULL)
  obs <- median(d$obs_min)
  nul <- vapply(seq_len(NPERM), function(i)
    median(vapply(d$m, function(mm) min(sample(pool, mm, replace=TRUE)), numeric(1))),
    numeric(1))
  tibble(a=d$a[1], b=d$b[1], loci=nrow(d), spec=round(median(pool),3),
         obs_min=round(obs,3), null_min=round(median(nul),3),
         null_lo=round(quantile(nul,.025),3), null_hi=round(quantile(nul,.975),3),
         ratio_naive=round(obs/median(pool),2),
         ratio_vs_null=round(obs/median(nul),3),
         p_two = 2*min(mean(nul >= obs), mean(nul <= obs)))
}))
out <- out %>% mutate(p_adj = round(p.adjust(p_two, "BH"), 4)) %>% arrange(ratio_vs_null)
print(as.data.frame(out), row.names=FALSE)

hr("2. verdict")
cat(sprintf("  pairs tested: %d\n", nrow(out)))
cat(sprintf("  obs within the null 95%% interval (SHARED)      : %d\n",
            sum(out$obs_min >= out$null_lo & out$obs_min <= out$null_hi)))
cat(sprintf("  obs ABOVE the null interval (INDEPENDENT)      : %d\n",
            sum(out$obs_min > out$null_hi)))
cat(sprintf("  obs BELOW the null interval (needs explaining) : %d\n",
            sum(out$obs_min < out$null_lo)))
cat(sprintf("\n  median naive ratio %.2f vs median ratio-to-null %.3f\n",
            median(out$ratio_naive), median(out$ratio_vs_null)))
cat("  The gap between those two columns is the estimator bias — the reason\n")
cat("  0.71-0.90 could not be read as evidence either way.\n")
write_csv(out, "trees/qc/shared_ab_null.csv")
