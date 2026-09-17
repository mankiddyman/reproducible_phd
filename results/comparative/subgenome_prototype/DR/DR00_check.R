#!/usr/bin/env Rscript
# DR00 check — what got built, what was lost, and DR01's scope.
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
LOC   <- read_tsv("DR/locus_meta.tsv", show_col_types=FALSE)
built <- sub("[.]fna$", "", list.files("DR/codon", pattern="[.]fna$"))
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",68), t, strrep("=",68)))

hr("1. build success by locus size")
f <- LOC %>% count(locus, name="tips") %>% mutate(built = locus %in% built)
print(as.data.frame(f %>% mutate(bin=cut(tips, c(0,2,3,5,10,50))) %>%
  group_by(bin) %>% summarise(loci=n(), built=sum(built),
                              failed=sum(!built), .groups="drop")), row.names=FALSE)
bad <- f %>% filter(!built, tips >= 3)
cat("\n  failures with >=3 tips (pal2nal, not size): n =", nrow(bad), "\n")
if (nrow(bad)) print(as.data.frame(LOC %>% distinct(locus, set, has_nep) %>%
  semi_join(bad, by="locus") %>% count(set, has_nep, name="failed")), row.names=FALSE)

hr("2. duplicated Dionaea genes (the many-to-many join)")
dup <- LOC %>% filter(genome=="Dionaea_muscipula") %>%
  count(tip, name="n_loci") %>% filter(n_loci > 1)
cat("  Dionaea genes appearing in >1 locus:", nrow(dup), "\n")
if (nrow(dup)) print(as.data.frame(dup %>% count(n_loci, name="genes")), row.names=FALSE)

hr("3. dS coverage")
k <- read_csv("DR/out/pairwise_ks.csv", show_col_types=FALSE)
cat("  rows:", nrow(k), " | loci with dS:", n_distinct(k$anchor), "\n")
cat("  usable (0<dS<3, codons>=100):",
    sum(k$dS>0 & k$dS<3 & k$codons>=100, na.rm=TRUE), "\n")
cat("  loci with >=1 usable pair:",
    n_distinct(k$anchor[k$dS>0 & k$dS<3 & k$codons>=100]), "\n")

hr("4. DR01 SCOPE — the decision point")
meta <- LOC %>% filter(locus %in% built)
cap <- meta %>% group_by(locus) %>%
  summarise(nep = any(genome=="Nepenthes_gracilis"),
            dio = sum(genome=="Dionaea_muscipula"), .groups="drop")
cat("  four-point capable (has Nepenthes):", sum(cap$nep), "\n")
cat("  delta capable (Dionaea == 2)      :", sum(cap$dio==2), "\n")
cat("  both                              :", sum(cap$nep & cap$dio==2), "\n\n")
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")
three <- meta %>% filter(genome %in% DROS) %>% count(locus, genome, name="k") %>%
  filter(k == 3) %>% semi_join(cap %>% filter(nep), by="locus") %>%
  mutate(region = sub("-.*$", "", locus))
old <- c(Drosera_regia=139, Drosera_capensis=105, Drosera_paradoxa=31,
         Drosera_binata=24, Drosera_scorpioides=16)
cat("  loci with exactly 3 copies of one species + a Nepenthes outgroup:\n")
print(as.data.frame(three %>% count(genome, name="loci_now") %>%
  mutate(loci_before = unname(old[genome]),
         gain = round(loci_now/loci_before, 1))), row.names=FALSE)
cat("\n  by ancestral region (chi-square at 2 df wants ~25 loci):\n")
print(as.data.frame(three %>% count(genome, region) %>% group_by(genome) %>%
  summarise(regions=n(), ge15=sum(n>=15), ge25=sum(n>=25),
            median_loci=median(n), .groups="drop")), row.names=FALSE)
cat("\n  VERDICT: species with >=3 regions at n>=25 are testable for AAB vs ABC.\n")
cat("  The rest inherit their constitution via d2/d1 in DR03.\n")
