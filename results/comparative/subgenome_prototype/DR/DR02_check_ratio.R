#!/usr/bin/env Rscript
# Is the ~1.1 pooled A:B ratio evidence against a 2:1 constitution?
# NO, if regions are individually 2:1 but disagree on WHICH side is doubled.
# Pooled counts cannot see that; per-region composition can.
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE)

cat("=== composition per (species, region) ===\n")
comp <- SS %>% filter(label %in% c("A","B")) %>%
  count(genome, region, label) %>%
  pivot_wider(names_from=label, values_from=n, values_fill=0) %>%
  mutate(tot=A+B, patt=paste0(A,"A:",B,"B"),
         is21=(pmax(A,B)==2 & pmin(A,B)==1),
         doubled=ifelse(A>B, "A", ifelse(B>A, "B", "tie")))
print(as.data.frame(comp %>% filter(tot==3) %>% count(genome, patt) %>%
  pivot_wider(names_from=patt, values_from=n, values_fill=0)), row.names=FALSE)

cat("\n=== among 3-segment regions ===\n")
t3 <- comp %>% filter(tot==3)
cat(sprintf("  regions with exactly 3 called segments: %d\n", nrow(t3)))
cat(sprintf("  of those, 2:1 in either direction: %d (%.0f%%)\n",
            sum(t3$is21), 100*mean(t3$is21)))
cat(sprintf("  3:0 or 0:3 (all one side): %d\n", sum(!t3$is21)))
cat("\n  A binomial null at the observed marginal A-rate would give ~44%% at\n")
cat("  exactly 2:1 (script 21). Well above that = the under-dispersion signal.\n")

cat("\n=== which side is doubled, per region? ===\n")
print(as.data.frame(t3 %>% count(genome, doubled) %>%
  pivot_wider(names_from=doubled, values_from=n, values_fill=0)), row.names=FALSE)
dA <- sum(t3$doubled=="A"); dB <- sum(t3$doubled=="B")
cat(sprintf("\n  A doubled in %d regions, B doubled in %d (%.0f%% favour A)\n",
            dA, dB, 100*dA/(dA+dB)))
cat(sprintf("  binomial p vs 50/50: %.4f\n", binom.test(dA, dA+dB, 0.5)$p.value))
cat("  script 36b attacks 6-7 gave an A share of 0.708 against a 0.500 null.\n")

cat("\n=== predicted pooled ratio from this orientation mix ===\n")
f <- dA/(dA+dB)
cat(sprintf("  if every region is 2:1 and %.2f favour A:\n", f))
cat(sprintf("    predicted pooled A:B = %.2f\n", (2*f+1*(1-f))/(1*f+2*(1-f))))
obs <- SS %>% filter(label %in% c("A","B")) %>% count(label) %>%
       pivot_wider(names_from=label, values_from=n)
cat(sprintf("    observed pooled A:B  = %.2f\n", obs$A/obs$B))
cat("\n  if these match, the pooled ratio is NOT evidence against 2:1 --\n")
cat("  it is 2:1 within regions plus inconsistent orientation between them.\n")

cat("\n=== do regions agree ACROSS species? ===\n")
cat("  Under shared ancestry, the same region should double the same side in\n")
cat("  every species. Disagreement means orientation is region-specific noise.\n\n")
print(as.data.frame(t3 %>% select(genome, region, doubled) %>%
  pivot_wider(names_from=genome, values_from=doubled)), row.names=FALSE)
agree <- t3 %>% count(region, doubled) %>% group_by(region) %>%
  summarise(species=sum(n), top=max(n), .groups="drop")
cat(sprintf("\n  mean within-region agreement across species: %.2f (0.5 = chance)\n",
            mean(agree$top/agree$species)))

cat("\n\n=== DOES AB03's HOMOEOLOGOUS EXCHANGE EXPLAIN THE RESIDUAL? ===\n")
cat("  AB02 found real exchange tracts in Dionaea: chr5 (one interior, 61\n")
cat("  genes) and chr8 (two telomere-proximal, 37 and 25). Within those the A\n")
cat("  chromosome carries B ancestry, so every Drosera label there is inverted.\n")
cat("  Dionaea chr5 -> region chr6_dom | Dionaea chr8 -> region chr10_dom.\n\n")
ff <- "AB/out/AB03_flagged_genes.csv"
if (!file.exists(ff)) { cat("  AB03 output not found -- run AB03_carryforward.R\n") } else {
  flag <- read_csv(ff, show_col_types=FALSE)
  CAL <- tibble(dpair=paste0("chr",1:8),
                region=c("chr2_dom","chr3_dom","chr4_dom","chr5_dom",
                         "chr6_dom","chr7_dom","chr9_dom","chr10_dom"))
  ft <- flag %>% filter(in_flip_region) %>% count(pair, name="flagged_loci") %>%
        left_join(CAL, by=c("pair"="dpair"))
  cat("  Dionaea exchange tracts, in Nepenthes region coordinates:\n")
  print(as.data.frame(ft), row.names=FALSE)
  cat("\n  do those regions still disagree across species after the fix?\n")
  chk <- t3 %>% filter(region %in% ft$region) %>%
         select(genome, region, doubled)
  if (nrow(chk)) print(as.data.frame(chk %>%
    pivot_wider(names_from=genome, values_from=doubled)), row.names=FALSE)
  cat("\n  A region still discordant AND carrying an exchange tract is\n")
  cat("  explained. Discordant with NO tract needs another explanation.\n")
}
