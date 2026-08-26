#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

h <- read_csv("synteny_ortho_hits.csv", show_col_types = FALSE)
k <- h %>% count(nep_gene, sp, name = "k") %>%
  pivot_wider(names_from = sp, values_from = k, values_fill = 0)

full <- k %>% filter(Dionaea_muscipula == 2, Drosera_regia == 3, Drosera_binata == 3)
cat(sprintf("FULL set (Dio=2, regia=3, binata=3): %d anchors\n", nrow(full)))

# near-complete fallback, for power
near <- k %>% filter(Dionaea_muscipula == 2,
                     Drosera_regia >= 2, Drosera_regia <= 3,
                     Drosera_binata >= 2, Drosera_binata <= 3) %>%
  anti_join(select(full, nep_gene), by = "nep_gene")
cat(sprintf("NEAR set (Dio=2, regia 2-3, binata 2-3): %d anchors\n", nrow(near)))

sel <- h %>% filter(sp %in% c("Dionaea_muscipula","Drosera_regia","Drosera_binata"),
                    nep_gene %in% c(full$nep_gene, near$nep_gene))

tipm <- bind_rows(
  sel %>% distinct(nep_gene, nep_chr) %>%
    transmute(nep_gene, tip = paste0("NEP@", nep_gene), gene = nep_gene,
              genome = "Nepenthes_gracilis", chr = nep_chr),
  sel %>% transmute(nep_gene, tip = paste0(sp, "@", sp_gene), gene = sp_gene,
                    genome = sp, chr = sp_chr)) %>%
  left_join(distinct(sel, nep_gene, nep_chr), by = "nep_gene") %>%
  mutate(set = ifelse(nep_gene %in% full$nep_gene, "full", "near"))

write_tsv(tipm, "full14/tip_meta.tsv")
for (a in unique(tipm$nep_gene))
  writeLines(tipm$tip[tipm$nep_gene == a], file.path("full14/ids", paste0(a, ".ids")))
writeLines(unique(filter(tipm, set == "full")$nep_gene), "full14/list_full.txt")
writeLines(unique(filter(tipm, set == "near")$nep_gene), "full14/list_near.txt")

cat("\n=== FULL anchors, per ancestral region ===\n")
print(as.data.frame(count(distinct(filter(tipm, set=="full"), nep_gene, nep_chr), nep_chr)))
