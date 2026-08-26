#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Cross-pair coherence screen.
# For HOGs with exactly one gene on each member of a Dionaea pair plus >=1 outgroup
# gene, ask which Dionaea copy the outgroup gene scores higher against.
# Drosera may be nested between the progenitors -> real topological signal.
# Nepenthes is outside both -> MUST vote at chance. That is the control.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
res <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE)
part <- res %>% filter(anchor == "Nepenthes_gracilis") %>%
  transmute(pair_lab, X = winner, Y = ifelse(winner == chrA, chrB, chrA))
cat("=== partition under test ===\n"); print(as.data.frame(part))

b <- bed %>% mutate(isRep = as.logical(isArrayRep)) %>% filter(is.na(isRep) | isRep)

xy <- bind_rows(transmute(part, chr = X, side = "X", pair_lab),
                transmute(part, chr = Y, side = "Y", pair_lab))
dio <- b %>% filter(genome == "Dionaea_muscipula") %>%
  inner_join(xy, by = "chr") %>% select(globHOG, dio_gene = id, side, pair_lab)

# HOGs with exactly one X copy and one Y copy from the SAME pair
clean <- dio %>% group_by(globHOG) %>%
  filter(n() == 2, n_distinct(side) == 2, n_distinct(pair_lab) == 1) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(globHOG, pair_lab), names_from = side, values_from = dio_gene)
cat(sprintf("\nclean 1X:1Y HOGs: %d\n", nrow(clean)))

read_of <- function(sp) {
  f <- file.path(GSD, sprintf("%s__v__Dionaea_muscipula.tsv", sp))
  if (!file.exists(f)) { message("missing: ", f); return(NULL) }
  x <- read.table(f, sep = "\t", header = FALSE, quote = "", comment.char = "",
                  stringsAsFactors = FALSE)
  # OrthoFinder BLAST-format tsv: q, s, pident, len, mm, go, qs, qe, ss, se, e, bits
  data.frame(q = x[[1]], s = x[[2]], bits = as.numeric(x[[ncol(x)]]),
             stringsAsFactors = FALSE)
}

vote <- function(sp) {
  h <- read_of(sp); if (is.null(h)) return(NULL)
  hx <- h %>% inner_join(select(clean, globHOG, pair_lab, X), by = c("s" = "X")) %>%
    group_by(q, globHOG, pair_lab) %>% summarise(bX = max(bits), .groups = "drop")
  hy <- h %>% inner_join(select(clean, globHOG, pair_lab, Y), by = c("s" = "Y")) %>%
    group_by(q, globHOG, pair_lab) %>% summarise(bY = max(bits), .groups = "drop")
  inner_join(hx, hy, by = c("q", "globHOG", "pair_lab")) %>%
    mutate(d = bX - bY, rel = d / pmax(bX, bY), sp = sp) %>%
    filter(rel != 0)
}

sps <- c(grep("^Drosera", unique(b$genome), value = TRUE), "Nepenthes_gracilis")
v <- bind_rows(lapply(sps, vote))
if (!nrow(v)) stop("no votes — check OrthoFinder tsv column layout")
write_csv(v, "topology_votes.csv")

cat("\n=== per species x pair: fraction of genes voting X ===\n")
pp <- v %>% group_by(sp, pair_lab) %>%
  summarise(n = n(), frac_X = mean(d > 0),
            med_rel = median(rel),
            p = binom.test(sum(d > 0), n(), 0.5)$p.value, .groups = "drop")
print(as.data.frame(pp %>% select(sp, pair_lab, n, frac_X, p) %>%
                      pivot_wider(names_from = pair_lab, values_from = c(frac_X, p)) %>%
                      select(sp, starts_with("frac_X"))), digits = 3)

cat("\n=== pooled per species ===\n")
sm <- v %>% group_by(sp) %>%
  summarise(n = n(), frac_X = mean(d > 0), med_rel = median(rel),
            p_pool = binom.test(sum(d > 0), n(), 0.5)$p.value, .groups = "drop") %>%
  left_join(pp %>% group_by(sp) %>%
              summarise(pairs_X = sum(frac_X > 0.5),
                        p_sign = binom.test(max(sum(frac_X > 0.5), 8 - sum(frac_X > 0.5)),
                                            8, 0.5)$p.value, .groups = "drop"), by = "sp")
print(as.data.frame(sm), digits = 3)

cat("\nREAD: Nepenthes frac_X must be ~0.5 (outside both progenitors).\n")
cat("      Drosera consistently off 0.5 in the SAME direction across pairs\n")
cat("      = cross-pair coherence = X and Y are ancestral genomes.\n")
cat("      Drosera ~ Nepenthes = no topological signal; WGD likely Dionaea-specific.\n")
cat("      Drosera off 0.5 but direction varies BY PAIR = chromosome-autonomous, no partition.\n")
