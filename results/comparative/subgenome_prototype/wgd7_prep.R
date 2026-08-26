#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Only two conditions: Nepenthes anchor present, Dionaea exactly 2 copies.
# NO cap on Drosera copy number -- capensis is a collapsed dodecaploid and can
# legitimately carry 6+ dispersed copies. Tandem duplicates were already removed
# upstream by the GENESPACE arrayID filter.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

h <- read_csv("synteny_ortho_hits.csv", show_col_types = FALSE)
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")
h <- filter(h, sp %in% c("Dionaea_muscipula", DROS))

k <- h %>% count(nep_gene, sp, name = "k") %>%
  pivot_wider(names_from = sp, values_from = k, values_fill = 0)
keep <- k %>% filter(Dionaea_muscipula == 2) %>% pull(nep_gene)
sel <- filter(h, nep_gene %in% keep)

cat(sprintf("anchors with Dio==2: %d\n", length(keep)))
cat("\ncopies per anchor, by species (no cap applied):\n")
print(as.data.frame(sel %>% count(nep_gene, sp, name="k") %>% group_by(sp) %>%
  summarise(anchors = n(), median = median(k), mean = round(mean(k),2),
            p90 = quantile(k, .9), max = max(k), .groups="drop")))
cat("\ntips per tree distribution:\n")
tt <- sel %>% count(nep_gene, name = "n") %>% mutate(n = n + 1)
print(summary(tt$n)); cat(sprintf("trees with >30 tips: %d | >50: %d\n",
                                  sum(tt$n > 30), sum(tt$n > 50)))

tips <- bind_rows(
  sel %>% distinct(nep_gene, nep_chr) %>%
    transmute(nep_gene, tip = paste0("Nepenthes_gracilis@", nep_gene),
              gene = nep_gene, genome = "Nepenthes_gracilis", chr = nep_chr),
  sel %>% transmute(nep_gene, tip = paste0(sp, "@", sp_gene),
                    gene = sp_gene, genome = sp, chr = sp_chr))
dir.create("wgd7/ids", recursive = TRUE, showWarnings = FALSE)
dir.create("wgd7/fa", showWarnings = FALSE); dir.create("wgd7/aln", showWarnings = FALSE)
dir.create("wgd7/tre", showWarnings = FALSE)
write_tsv(tips, "wgd7/tip_meta.tsv")
for (a in unique(tips$nep_gene))
  writeLines(tips$tip[tips$nep_gene == a], file.path("wgd7/ids", paste0(a, ".ids")))
writeLines(unique(tips$nep_gene), "wgd7/anchorlist.txt")
cat(sprintf("\nwrote %d anchors | %d tips | median %.0f tips/tree\n",
            n_distinct(tips$nep_gene), nrow(tips), median(table(tips$nep_gene))))
print(as.data.frame(count(tips, genome)))
