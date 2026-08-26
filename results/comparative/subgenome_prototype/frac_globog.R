#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Same fractionation test, globOG (flat orthogroups) instead of globHOG.
# globOG does not apply the hierarchical duplication-node cut that splits
# rate-accelerated homeologs, so it should recover copies HOGs lose.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
b <- bed %>% mutate(isRep = as.logical(isArrayRep)) %>% filter(is.na(isRep) | isRep)
ref <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE) %>%
  filter(anchor == "Nepenthes_gracilis") %>% select(pair_lab, chrA, chrB,
                                                    hog_frac = frac_A, hog_winner = winner)

run <- function(col, tag) {
  bb <- b %>% rename(grp = all_of(col))
  dio <- bb %>% filter(genome == "Dionaea_muscipula") %>%
    mutate(pl = sub("_sg[0-9]+_s[0-9]+$", "", chr)) %>% select(grp, chr, pl)
  st <- dio %>% group_by(grp) %>%
    summarise(n_dio = n(), pl = pl[1], one_chr = ifelse(n() == 1, chr[1], NA), .groups = "drop")
  anc <- bb %>% filter(genome != "Dionaea_muscipula") %>% distinct(grp) %>% pull(grp)
  lost <- st %>% filter(n_dio == 1, grp %in% anc, !is.na(one_chr))
  ref %>% rowwise() %>%
    mutate(kA = sum(lost$one_chr == chrA), kB = sum(lost$one_chr == chrB)) %>% ungroup() %>%
    mutate(total = kA + kB, frac_A = kA/total,
           p = mapply(function(a,t) if (t>=5) binom.test(a,t,0.5)$p.value else NA, kA, total),
           p_adj = p.adjust(p, "BH"),
           winner = ifelse(kA > kB, chrA, chrB), src = tag) %>%
    select(pair_lab, src, kA, kB, total, frac_A, p_adj, winner)
}

hg <- run("globHOG", "globHOG"); og <- run("globOG", "globOG")
cmp <- inner_join(select(hg, pair_lab, hog = frac_A, n_hog = total, w_hog = winner),
                  select(og, pair_lab, og = frac_A, n_og = total, w_og = winner),
                  by = "pair_lab") %>%
  left_join(select(ref, pair_lab, hog_winner), by = "pair_lab") %>%
  mutate(agree = w_og == w_hog,
         dev_hog = abs(hog - .5), dev_og = abs(og - .5))

cat("=== fractionation: globHOG vs globOG ===\n")
print(as.data.frame(cmp %>% select(pair_lab, n_hog, hog, n_og, og, w_og, agree)), digits = 3)
cat(sprintf("\ndirection agrees %d/8 | mean |dev| HOG %.3f vs OG %.3f\n",
            sum(cmp$agree), mean(cmp$dev_hog), mean(cmp$dev_og)))
cat(sprintf("usable groups: HOG %d, OG %d\n", sum(cmp$n_hog), sum(cmp$n_og)))
write_csv(bind_rows(hg, og), "fractionation_og_vs_hog.csv")
