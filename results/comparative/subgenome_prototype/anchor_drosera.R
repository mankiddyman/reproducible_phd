#!/usr/bin/env Rscript
# PIN LIBPATH: genespace module sets R_LIBS_USER (rlang 1.1.3, ggplot2 3.5.1)
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))

# Synteny-free fractionation replication across 6 anchors.
# Anchor certifies ancestral presence only. Dionaea's own chromosome gives the pair.
# Immune to syntenic-block detection artifacts by construction.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)

cat("=== genomes in GENESPACE run ===\n")
print(bed %>% count(genome, name = "n_genes") %>% as.data.frame())

b <- bed %>% mutate(isRep = as.logical(isArrayRep)) %>% filter(is.na(isRep) | isRep)

pairs <- b %>% filter(genome == "Dionaea_muscipula") %>%
  mutate(pair_lab = sub("_sg[0-9]+_s[0-9]+$", "", chr)) %>%
  distinct(pair_lab, chr) %>% arrange(pair_lab, chr) %>%
  group_by(pair_lab) %>%
  summarise(chrA = chr[1], chrB = chr[2], n_chr = n(), .groups = "drop")
stopifnot(all(pairs$n_chr == 2))

dio <- b %>% filter(genome == "Dionaea_muscipula") %>%
  mutate(pair_lab = sub("_sg[0-9]+_s[0-9]+$", "", chr)) %>%
  select(globHOG, dio_chr = chr, pair_lab)

state <- dio %>% group_by(globHOG) %>%
  summarise(n_dio = n(), n_pair = n_distinct(pair_lab), n_chr = n_distinct(dio_chr),
            pair_lab = pair_lab[1], one_chr = ifelse(n() == 1, dio_chr[1], NA_character_),
            .groups = "drop")

lost1 <- state %>% filter(n_dio == 1)                              # one homeolog gone
ret2  <- state %>% filter(n_dio == 2, n_pair == 1, n_chr == 2)     # both kept, clean pair

cat(sprintf("\nDionaea HOGs: 1 copy %d | 2 copies clean-pair %d | other %d\n",
            nrow(lost1), nrow(ret2), nrow(state) - nrow(lost1) - nrow(ret2)))

anchors <- setdiff(unique(b$genome), "Dionaea_muscipula")

run <- function(sp) {
  hg <- b %>% filter(genome == sp) %>% distinct(globHOG) %>% pull(globHOG)
  L <- lost1 %>% filter(globHOG %in% hg, !is.na(pair_lab))
  R <- ret2  %>% filter(globHOG %in% hg)
  if (!nrow(L)) return(NULL)
  out <- pairs %>% select(pair_lab, chrA, chrB) %>%
    rowwise() %>%
    mutate(kA = sum(L$one_chr == chrA), kB = sum(L$one_chr == chrB)) %>%
    ungroup() %>%
    mutate(total = kA + kB, frac_A = kA / total,
           p = mapply(function(a, t) if (t >= 5) binom.test(a, t, 0.5)$p.value else NA_real_,
                      kA, total),
           p_adj = p.adjust(p, "BH"),
           winner = ifelse(kA > kB, chrA, chrB),
           anchor = sp,
           n_ret = sapply(pair_lab, function(p) sum(R$pair_lab == p)))
  out
}

res <- bind_rows(lapply(anchors, run))
write_csv(res, "anchor_drosera_fractionation.csv")

for (sp in anchors) {
  r <- filter(res, anchor == sp); if (!nrow(r)) next
  cat(sprintf("\n=== %s ===\n", sp))
  print(as.data.frame(select(r, pair_lab, kA, kB, total, frac_A, p_adj, winner, n_ret)),
        digits = 3)
}

## ---- sign-vector agreement matrix -------------------------------------------
wide <- res %>% select(anchor, pair_lab, winner) %>%
  pivot_wider(names_from = anchor, values_from = winner) %>% arrange(pair_lab)
cat("\n=== winner per pair, per anchor ===\n"); print(as.data.frame(wide))

ref <- read_csv("fractionation_by_chrpair.csv", show_col_types = FALSE) %>%
  select(pair_lab = exp_pair, ref_winner = retained_more)
agr <- wide %>% left_join(ref, by = "pair_lab")
cat("\n=== agreement with Nepenthes-synteny result ===\n")
for (sp in anchors) {
  if (!sp %in% names(agr)) next
  k <- sum(agr[[sp]] == agr$ref_winner, na.rm = TRUE); n <- sum(!is.na(agr[[sp]]))
  cat(sprintf("%-24s %d/%d  sign-test p=%.4f\n", sp, k, n,
              binom.test(max(k, n - k), n, 0.5)$p.value))
}

## ---- majority vote across the 5 Drosera anchors -----------------------------
dros <- grep("^Drosera", anchors, value = TRUE)
if (length(dros) >= 2) {
  mv <- res %>% filter(anchor %in% dros) %>% count(pair_lab, winner) %>%
    group_by(pair_lab) %>% slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>% rename(drosera_vote = winner, n_anchors = n) %>%
    left_join(ref, by = "pair_lab") %>%
    mutate(agree = drosera_vote == ref_winner)
  cat("\n=== Drosera majority vote vs Nepenthes ===\n")
  print(as.data.frame(mv))
  k <- sum(mv$agree, na.rm = TRUE); n <- nrow(mv)
  cat(sprintf("agree %d/%d | sign-test p=%.4f | floor with one label fixed %.4f\n",
              k, n, binom.test(max(k, n - k), n, 0.5)$p.value, 0.5^(n - 1)))
  write_csv(mv, "drosera_majority_vote.csv")
}

## ---- how independent is the evidence? ---------------------------------------
cat("\n=== HOG-set overlap between anchors (Jaccard on loss events) ===\n")
sets <- lapply(anchors, function(sp) {
  hg <- b %>% filter(genome == sp) %>% distinct(globHOG) %>% pull(globHOG)
  lost1 %>% filter(globHOG %in% hg, !is.na(pair_lab)) %>% pull(globHOG)
})
names(sets) <- anchors
J <- outer(seq_along(sets), seq_along(sets), Vectorize(function(i, j)
  length(intersect(sets[[i]], sets[[j]])) / length(union(sets[[i]], sets[[j]]))))
dimnames(J) <- list(anchors, anchors)
print(round(J, 3))
cat("\nhigh Jaccard vs Nepenthes = same genes = weak replication\n")
