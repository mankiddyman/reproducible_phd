#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# 1. Gene-level effective ploidy: how many chromosomes carry a syntenic ortholog?
# 2. Are the k chromosomes per Nepenthes region genuine homeologs or minor hits?
# 3. Are chromosomes subgenome-coherent, or mosaics needing block-level assignment?
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

blk <- read_csv(file.path(GSD, "syntenicBlock_coordinates.csv"), show_col_types = FALSE)
bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)

nd <- blk %>% filter(genome1 == "Nepenthes_gracilis" | genome2 == "Nepenthes_gracilis") %>%
  mutate(n1 = genome1 == "Nepenthes_gracilis",
         sp = ifelse(n1, genome2, genome1),
         nep_chr = ifelse(n1, chr1, chr2),
         sp_chr  = ifelse(n1, chr2, chr1),
         nep_s = ifelse(n1, pmin(startBp1, endBp1), pmin(startBp2, endBp2)),
         nep_e = ifelse(n1, pmax(startBp1, endBp1), pmax(startBp2, endBp2)),
         nhit  = ifelse(n1, nHits1, nHits2)) %>%
  filter(sp != "Nepenthes_gracilis", grepl("_dom$", nep_chr), !is.na(nep_s))

merge_iv <- function(s, e) { o <- order(s); s <- s[o]; e <- e[o]
  cs <- s[1]; ce <- e[1]; tot <- 0
  if (length(s) > 1) for (i in 2:length(s)) {
    if (s[i] <= ce) ce <- max(ce, e[i]) else { tot <- tot + ce - cs; cs <- s[i]; ce <- e[i] } }
  tot + ce - cs }

cov <- nd %>% group_by(sp, nep_chr, sp_chr) %>%
  summarise(n_blk = n(), hits = sum(nhit, na.rm = TRUE),
            cov_Mb = merge_iv(nep_s, nep_e)/1e6, .groups = "drop") %>%
  group_by(sp, nep_chr) %>% mutate(share = cov_Mb/sum(cov_Mb)) %>% ungroup()
write_csv(cov, "block_coverage_matrix.csv")

cat("=== major vs minor tracts (major = >=10% of the region's covered bp) ===\n")
print(as.data.frame(cov %>% mutate(major = share >= 0.10) %>%
  group_by(sp, nep_chr) %>%
  summarise(n_major = sum(major), n_minor = sum(!major), .groups = "drop") %>%
  group_by(sp) %>% summarise(median_major = median(n_major),
                             range_major = paste(range(n_major), collapse = "-"),
                             median_minor = median(n_minor), .groups = "drop")))
cat("median_major = effective syntenic ploidy after dropping fragments\n")

cat("\n=== chromosome coherence: how many ancestral regions per chromosome? ===\n")
coh <- cov %>% filter(share >= 0.10) %>% count(sp, sp_chr, name = "n_regions")
print(as.data.frame(coh %>% group_by(sp) %>%
  summarise(n_chr = n(), pct_single_region = round(100*mean(n_regions == 1)),
            median_regions = median(n_regions),
            max_regions = max(n_regions), .groups = "drop")))
cat("pct_single_region high => chromosome-level assignment OK\n")
cat("pct_single_region low  => chromosomes are mosaics, assign at BLOCK level\n")

## gene-level effective ploidy via syntenic hits
hits <- read.table(gzfile(file.path(GSD, "allInBlkOgHits.txt.gz")), header = TRUE,
                   sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
cat(sprintf("\nallInBlkOgHits columns: %s\n", paste(names(hits), collapse = ", ")))
cat(sprintf("rows: %d\n", nrow(hits)))
