#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# How many orthogroups can support a subgenome-aware gene tree,
# and what chromosome-level resolution do we have per species?
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
b <- bed %>% mutate(isRep = as.logical(isArrayRep)) %>% filter(is.na(isRep) | isRep)

cat("=== chromosome inventory per genome ===\n")
inv <- b %>% group_by(genome) %>%
  summarise(n_genes = n(), n_chr = n_distinct(chr),
            chrs = paste(sort(unique(chr))[1:min(6, n_distinct(chr))], collapse = " "),
            .groups = "drop")
print(as.data.frame(inv))

cat("\n=== copies per HOG, by species (median / mode) ===\n")
cn <- b %>% count(globHOG, genome, name = "k")
print(as.data.frame(cn %>% group_by(genome) %>%
  summarise(n_hog = n(), median_k = median(k), max_k = max(k),
            pct_1 = round(100*mean(k == 1)), pct_2 = round(100*mean(k == 2)),
            pct_3 = round(100*mean(k == 3)), pct_ge4 = round(100*mean(k >= 4)),
            .groups = "drop")))

wide <- cn %>% pivot_wider(names_from = genome, values_from = k, values_fill = 0)

need <- function(df, ...) { f <- rlang::quos(...); filter(df, !!!f) %>% nrow() }
cat("\n=== orthogroups available under different tip requirements ===\n")
combos <- tribble(
  ~label, ~n,
  "Nep>=1, Dio==2",
    sum(wide$Nepenthes_gracilis >= 1 & wide$Dionaea_muscipula == 2),
  "Nep>=1, Dio==2, regia>=2",
    sum(wide$Nepenthes_gracilis >= 1 & wide$Dionaea_muscipula == 2 & wide$Drosera_regia >= 2),
  "Nep>=1, Dio==2, binata>=2",
    sum(wide$Nepenthes_gracilis >= 1 & wide$Dionaea_muscipula == 2 & wide$Drosera_binata >= 2),
  "Nep>=1, Dio==2, regia>=2, binata>=2",
    sum(wide$Nepenthes_gracilis >= 1 & wide$Dionaea_muscipula == 2 &
        wide$Drosera_regia >= 2 & wide$Drosera_binata >= 2),
  "Nep>=1, Dio==2, regia==3, binata==3",
    sum(wide$Nepenthes_gracilis >= 1 & wide$Dionaea_muscipula == 2 &
        wide$Drosera_regia == 3 & wide$Drosera_binata == 3),
  "all 7 genomes present, <=4 copies each",
    sum(rowSums(select(wide, -globHOG) >= 1) == (ncol(wide)-1) &
        rowSums(select(wide, -globHOG) > 4) == 0)
)
print(as.data.frame(combos))

## how well does chromosome partition the copies in the richest set?
rich <- wide %>%
  filter(Nepenthes_gracilis >= 1, Dionaea_muscipula == 2,
         Drosera_regia >= 2, Drosera_binata >= 2) %>% pull(globHOG)
cat(sprintf("\nrich set: %d HOGs\n", length(rich)))
for (sp in c("Drosera_regia","Drosera_binata")) {
  x <- b %>% filter(globHOG %in% rich, genome == sp) %>%
    group_by(globHOG) %>% summarise(k = n(), nchr = n_distinct(chr), .groups = "drop")
  cat(sprintf("%-18s copies on distinct chromosomes in %.1f%% of HOGs (median %d copies)\n",
              sp, 100*mean(x$nchr == x$k), median(x$k)))
}
write_csv(wide %>% filter(globHOG %in% rich), "tree_candidate_hogs.csv")
cat("\nWROTE: tree_candidate_hogs.csv\n")
