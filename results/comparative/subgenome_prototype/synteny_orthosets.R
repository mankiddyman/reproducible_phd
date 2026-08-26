#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Positional ortholog sets from syntenic gene hits (allInBlkOgHits), Nepenthes-anchored.
# Immune to the rate-acceleration problem that breaks tree-based orthogroup inference.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
map <- bed %>% select(ofID, gene = id, genome, chr, ord)

h <- read.table(gzfile(file.path(GSD, "allInBlkOgHits.txt.gz")), header = TRUE,
                sep = "\t", stringsAsFactors = FALSE)
h <- bind_rows(transmute(h, a = ofID1, b = ofID2),
               transmute(h, a = ofID2, b = ofID1)) %>% distinct()
cat(sprintf("directed syntenic hits: %d\n", nrow(h)))

h <- h %>% inner_join(map, by = c("a" = "ofID")) %>%
  rename(nep_gene = gene, nep_gen = genome, nep_chr = chr, nep_ord = ord) %>%
  inner_join(map, by = c("b" = "ofID")) %>%
  rename(sp_gene = gene, sp = genome, sp_chr = chr, sp_ord = ord) %>%
  filter(nep_gen == "Nepenthes_gracilis", grepl("_dom$", nep_chr),
         sp != "Nepenthes_gracilis")
cat(sprintf("Nepenthes-anchored hits: %d | anchor genes: %d\n",
            nrow(h), n_distinct(h$nep_gene)))
write_csv(h, "synteny_ortho_hits.csv")

k <- h %>% count(nep_gene, sp, name = "k")
cat("\n=== syntenic copies per Nepenthes anchor gene ===\n")
print(as.data.frame(k %>% group_by(sp) %>%
  summarise(anchors = n(), median_k = median(k), mean_k = round(mean(k), 2),
            pct_1 = round(100*mean(k == 1)), pct_2 = round(100*mean(k == 2)),
            pct_3 = round(100*mean(k == 3)), pct_ge4 = round(100*mean(k >= 4)),
            .groups = "drop")))

w <- k %>% pivot_wider(names_from = sp, values_from = k, values_fill = 0)
cat("\n=== anchors usable for each monophyly test ===\n")
tests <- tribble(~test, ~n,
  "Dionaea WGD: Dio==2 + any Drosera>=1",
    sum(w$Dionaea_muscipula == 2 & (w$Drosera_regia >= 1 | w$Drosera_binata >= 1)),
  "regia WGD: regia>=2 + Dio>=1",
    sum(w$Drosera_regia >= 2 & w$Dionaea_muscipula >= 1),
  "binata WGD: binata>=2 + Dio>=1",
    sum(w$Drosera_binata >= 2 & w$Dionaea_muscipula >= 1),
  "regia vs binata: both >=2",
    sum(w$Drosera_regia >= 2 & w$Drosera_binata >= 2),
  "core set: Dio==2, regia>=2, binata>=2",
    sum(w$Dionaea_muscipula == 2 & w$Drosera_regia >= 2 & w$Drosera_binata >= 2),
  "full: Dio==2, regia==3, binata==3",
    sum(w$Dionaea_muscipula == 2 & w$Drosera_regia == 3 & w$Drosera_binata == 3))
print(as.data.frame(tests))

core <- w %>% filter(Dionaea_muscipula == 2, Drosera_regia >= 2, Drosera_binata >= 2,
                     Drosera_regia <= 4, Drosera_binata <= 4)
write_csv(core, "orthoset_core.csv")
cat(sprintf("\ncore set written: %d anchors, median tips %.0f\n", nrow(core),
            median(1 + core$Dionaea_muscipula + core$Drosera_regia + core$Drosera_binata)))

## does the Dionaea X/Y partition hold at gene level in this set?
ref <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE) %>%
  filter(anchor == "Nepenthes_gracilis")
xy <- bind_rows(transmute(ref, chr = winner, side = "X"),
                transmute(ref, chr = ifelse(winner == chrA, chrB, chrA), side = "Y"))
chk <- h %>% filter(sp == "Dionaea_muscipula", nep_gene %in% core$nep_gene) %>%
  left_join(xy, by = c("sp_chr" = "chr")) %>%
  group_by(nep_gene) %>% summarise(sides = paste(sort(side), collapse = ""), .groups = "drop")
cat("\n=== Dionaea copies per anchor, by X/Y side ===\n")
print(as.data.frame(count(chk, sides)))
cat("XY = one from each subgenome (usable); XX or YY = both from one side\n")
