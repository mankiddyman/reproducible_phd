#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Does the retention asymmetry survive normalising by each chromosome's gene count?
# If frac_A tracks gene-count share exactly, the HOG analysis adds nothing beyond
# "chrX is annotated more densely" and the whole result rests on annotation quality.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
FAI <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/Dionaea_muscipula/assembly_final/external_collapsed/Dionaea_muscipula_chr.fa.fai"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
res <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE)
fai <- read.table(FAI, stringsAsFactors = FALSE)[,1:2]; names(fai) <- c("chr","len")

b <- bed %>% mutate(isRep = as.logical(isArrayRep)) %>% filter(is.na(isRep) | isRep)
dm <- b %>% filter(genome == "Dionaea_muscipula")

## per-chromosome annotation profile
prof <- dm %>% group_by(chr) %>%
  summarise(n_genes = n(),
            n_inHOG = sum(!is.na(globHOG) & globHOG != ""),
            frac_inHOG = n_inHOG / n_genes,
            med_pepLen = median(pepLen, na.rm = TRUE),
            frac_noAnchor = mean(as.logical(noAnchor), na.rm = TRUE),
            .groups = "drop") %>%
  left_join(fai, by = "chr") %>% mutate(dens = n_genes / (len/1e6))
cat("=== Dionaea per-chromosome annotation profile ===\n")
print(as.data.frame(prof), digits = 3)

## no-anchor baseline: single-copy Dionaea HOGs, no certification at all
dio <- dm %>% mutate(pair_lab = sub("_sg[0-9]+_s[0-9]+$", "", chr))
st <- dio %>% group_by(globHOG) %>%
  summarise(n_dio = n(), pl = pair_lab[1], one_chr = ifelse(n() == 1, chr[1], NA),
            .groups = "drop") %>% filter(n_dio == 1, !is.na(one_chr))

pairs <- dio %>% distinct(pair_lab, chr) %>% arrange(pair_lab, chr) %>%
  group_by(pair_lab) %>% summarise(chrA = chr[1], chrB = chr[2], .groups = "drop")

base <- pairs %>% rowwise() %>%
  mutate(kA = sum(st$one_chr == chrA), kB = sum(st$one_chr == chrB)) %>% ungroup() %>%
  mutate(frac_A_noanchor = kA/(kA+kB))

## expected frac from raw gene counts alone
cmp <- pairs %>%
  left_join(select(prof, chr, gA = n_genes, hA = n_inHOG), by = c("chrA"="chr")) %>%
  left_join(select(prof, chr, gB = n_genes, hB = n_inHOG), by = c("chrB"="chr")) %>%
  mutate(exp_from_genes = gA/(gA+gB), exp_from_HOGgenes = hA/(hA+hB)) %>%
  left_join(select(base, pair_lab, frac_A_noanchor), by = "pair_lab") %>%
  left_join(res %>% filter(anchor == "Nepenthes_gracilis") %>%
              select(pair_lab, frac_nep = frac_A), by = "pair_lab") %>%
  left_join(res %>% filter(anchor == "Drosera_scorpioides") %>%
              select(pair_lab, frac_scorp = frac_A), by = "pair_lab") %>%
  mutate(excess_nep = frac_nep - exp_from_genes)

cat("\n=== observed retention vs gene-count expectation ===\n")
print(as.data.frame(select(cmp, pair_lab, exp_from_genes, exp_from_HOGgenes,
                           frac_A_noanchor, frac_nep, frac_scorp, excess_nep)), digits = 3)
c1 <- with(cmp, cor.test(exp_from_genes, frac_nep, method = "spearman"))
cat(sprintf("\nSpearman(gene-count share, anchored retention) rho=%.3f p=%.4f\n",
            c1$estimate, c1$p.value))
cat(sprintf("mean |excess over gene-count expectation| = %.4f\n", mean(abs(cmp$excess_nep))))
cat(sprintf("no-anchor vs Nepenthes-anchored direction agrees %d/8\n",
            sum((cmp$frac_A_noanchor > .5) == (cmp$frac_nep > .5))))

## is the X set under-annotated on any axis?
part <- res %>% filter(anchor == "Nepenthes_gracilis") %>%
  transmute(pair_lab, X = winner, Y = ifelse(winner == chrA, chrB, chrA))
qc <- part %>%
  left_join(select(prof, chr, hX = frac_inHOG, pX = med_pepLen, aX = frac_noAnchor),
            by = c("X"="chr")) %>%
  left_join(select(prof, chr, hY = frac_inHOG, pY = med_pepLen, aY = frac_noAnchor),
            by = c("Y"="chr"))
cat("\n=== annotation-quality axes, X vs Y ===\n"); print(as.data.frame(qc), digits = 3)
for (v in list(c("hX","hY","frac_inHOG"), c("pX","pY","med_pepLen"),
               c("aX","aY","frac_noAnchor"))) {
  x <- qc[[v[1]]]; y <- qc[[v[2]]]
  cat(sprintf("%-14s X>Y in %d/8 | paired Wilcoxon p=%.3f\n", v[3], sum(x > y, na.rm=TRUE),
              suppressWarnings(wilcox.test(x, y, paired = TRUE)$p.value)))
}
write_csv(cmp, "annotation_reduction.csv")
