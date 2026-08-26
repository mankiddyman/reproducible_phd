#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Same design as wgd_prep.R but with ALL FIVE Drosera species, not just regia + binata.
# Dionaea must have exactly 2 copies (the phasing unit); Drosera capped at 4 copies each.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
CAP <- 4; MAXTIPS <- 18

h <- read_csv("synteny_ortho_hits.csv", show_col_types = FALSE)
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")
h <- filter(h, sp %in% c("Dionaea_muscipula", DROS))

k <- h %>% count(nep_gene, sp, name = "k")
ok <- k %>% filter(k <= CAP) %>% group_by(nep_gene) %>%
  filter(any(sp == "Dionaea_muscipula" & k == 2), sum(k) + 1 <= MAXTIPS) %>% ungroup()
sel <- h %>% semi_join(distinct(ok, nep_gene, sp), by = c("nep_gene","sp")) %>%
  group_by(nep_gene) %>% filter(any(sp == "Dionaea_muscipula")) %>% ungroup()

tips <- bind_rows(
  sel %>% distinct(nep_gene, nep_chr) %>%
    transmute(nep_gene, tip = paste0("Nepenthes_gracilis@", nep_gene),
              gene = nep_gene, genome = "Nepenthes_gracilis", chr = nep_chr),
  sel %>% transmute(nep_gene, tip = paste0(sp, "@", sp_gene),
                    gene = sp_gene, genome = sp, chr = sp_chr))
write_tsv(tips, "wgd6/tip_meta.tsv")
for (a in unique(tips$nep_gene))
  writeLines(tips$tip[tips$nep_gene == a], file.path("wgd6/ids", paste0(a, ".ids")))
writeLines(unique(tips$nep_gene), "wgd6/anchorlist.txt")

cat(sprintf("anchors: %d | tips: %d | median tips/tree %.0f\n",
            n_distinct(tips$nep_gene), nrow(tips), median(table(tips$nep_gene))))
cat("\ntips per genome:\n"); print(as.data.frame(count(tips, genome)))
cat("\n--- for comparison, the 3-species run ---\n")
old <- read_tsv("wgd/tip_meta.tsv", show_col_types = FALSE)
cat(sprintf("anchors: %d | tips: %d\n", n_distinct(old$nep_gene), nrow(old)))
