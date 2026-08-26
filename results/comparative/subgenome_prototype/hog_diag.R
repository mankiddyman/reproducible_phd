#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# 1. Are homeologs being split across HOGs? Compare globOG vs globHOG copy counts.
# 2. Synteny is ground truth: how many chromosomes of each species hit one Nepenthes chr?
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
b <- bed %>% mutate(isRep = as.logical(isArrayRep)) %>% filter(is.na(isRep) | isRep)

cat("=== copies per group: globOG (flat) vs globHOG (hierarchical) ===\n")
f <- function(col) b %>% count(.data[[col]], genome, name = "k") %>%
  group_by(genome) %>% summarise(groups = n(), mean_k = mean(k),
                                 pct_1 = round(100*mean(k == 1)),
                                 pct_ge2 = round(100*mean(k >= 2)),
                                 pct_ge3 = round(100*mean(k >= 3)), .groups = "drop")
og <- f("globOG") %>% mutate(src = "globOG"); hg <- f("globHOG") %>% mutate(src = "globHOG")
print(as.data.frame(bind_rows(og, hg) %>% arrange(genome, src)), digits = 3)
cat("\nif globOG shows the copies and globHOG does not, HOG inference is splitting homeologs\n")

## synteny: chromosomes per Nepenthes chromosome
blk <- read_csv(file.path(GSD, "syntenicBlock_coordinates.csv"), show_col_types = FALSE)
nd <- blk %>%
  filter(genome1 == "Nepenthes_gracilis" | genome2 == "Nepenthes_gracilis") %>%
  mutate(n1 = genome1 == "Nepenthes_gracilis",
         sp = ifelse(n1, genome2, genome1),
         nep_chr = ifelse(n1, chr1, chr2),
         sp_chr = ifelse(n1, chr2, chr1)) %>%
  filter(sp != "Nepenthes_gracilis", grepl("_dom$", nep_chr))

cat("\n=== distinct chromosomes hitting each Nepenthes chromosome ===\n")
print(as.data.frame(nd %>% group_by(sp, nep_chr) %>%
  summarise(n_chr = n_distinct(sp_chr), n_blk = n(), .groups = "drop") %>%
  group_by(sp) %>% summarise(median_chr_per_nep = median(n_chr),
                             range = paste(range(n_chr), collapse = "-"),
                             total_blocks = sum(n_blk), .groups = "drop")))
cat("\nthis is the effective syntenic ploidy GENESPACE actually resolved\n")

cat("\n=== per-species detail: chromosomes grouped by Nepenthes affinity ===\n")
det <- nd %>% count(sp, nep_chr, sp_chr, name = "n_blk") %>%
  group_by(sp, nep_chr) %>% arrange(desc(n_blk), .by_group = TRUE) %>%
  summarise(chrs = paste(sp_chr, collapse = ", "), .groups = "drop")
for (s in unique(det$sp)) {
  cat(sprintf("\n-- %s --\n", s))
  print(as.data.frame(filter(det, sp == s) %>% select(nep_chr, chrs)), right = FALSE)
}
write_csv(det, "nepenthes_affinity_groups.csv")
