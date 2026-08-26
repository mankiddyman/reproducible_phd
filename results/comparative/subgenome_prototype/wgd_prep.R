#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# One alignment per Nepenthes anchor; all three monophyly tests read off the same trees.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
CAP <- 4

h <- read_csv("synteny_ortho_hits.csv", show_col_types = FALSE) %>%
  filter(sp %in% c("Dionaea_muscipula", "Drosera_regia", "Drosera_binata"))

ok <- h %>% count(nep_gene, sp, name = "k") %>% filter(k <= CAP) %>%
  group_by(nep_gene) %>%
  filter(n_distinct(sp) >= 2, any(k >= 2), sum(k) >= 3, sum(k) <= 10) %>%
  ungroup()

sel <- h %>% semi_join(distinct(ok, nep_gene, sp), by = c("nep_gene", "sp"))

tips <- bind_rows(
  sel %>% distinct(nep_gene) %>%
    mutate(tip = paste0("Nepenthes_gracilis@", nep_gene), gene = nep_gene,
           genome = "Nepenthes_gracilis", chr = NA_character_),
  sel %>% transmute(nep_gene, tip = paste0(sp, "@", sp_gene), gene = sp_gene,
                    genome = sp, chr = sp_chr))

write_tsv(tips, "wgd/tip_meta.tsv")
for (a in unique(tips$nep_gene))
  writeLines(tips$tip[tips$nep_gene == a], file.path("wgd/ids", paste0(a, ".ids")))
writeLines(unique(tips$nep_gene), "wgd/anchorlist.txt")

cat(sprintf("anchors: %d | tips: %d | median tips/tree %.0f\n",
            n_distinct(tips$nep_gene), nrow(tips), median(table(tips$nep_gene))))
cat("\ntips per genome:\n"); print(as.data.frame(count(tips, genome)))
