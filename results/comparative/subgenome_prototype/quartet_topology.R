#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Restrict each tree to a clean 5-taxon quartet: NEP + DioX + DioY + 2 copies of one
# other species on DISTINCT chromosomes. Two outcomes:
#   Dio monophyletic -> independent WGDs
#   interleaved      -> shared WGD, AND the pairing phases both species at once
# NULL for monophyly is tree-shape dependent (1/3 balanced, 1/6 pectinate) -> permute.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(ape) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

meta <- read_tsv("wgd/tip_meta.tsv", show_col_types = FALSE)
hits <- read_csv("synteny_ortho_hits.csv", show_col_types = FALSE) %>% distinct(nep_gene, nep_chr)
ref  <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE) %>%
  filter(anchor == "Nepenthes_gracilis")
xy <- bind_rows(transmute(ref, chr = winner, side = "X"),
                transmute(ref, chr = ifelse(winner == chrA, chrB, chrA), side = "Y"))
lab <- meta %>% left_join(xy, by = "chr")

fs <- list.files("wgd/tre", pattern = "\\.tre$", full.names = TRUE)
cat(sprintf("trees: %d\n", length(fs)))

one <- function(f) {
  tr <- tryCatch(read.tree(f), error = function(e) NULL)
  if (is.null(tr) || Ntip(tr) < 5) return(NULL)
  m <- lab[match(tr$tip.label, lab$tip), ]
  if (any(is.na(m$genome))) return(NULL)
  nep <- tr$tip.label[m$genome == "Nepenthes_gracilis"]
  if (length(nep) != 1) return(NULL)
  di <- which(m$genome == "Dionaea_muscipula")
  if (length(di) != 2 || !setequal(m$side[di], c("X", "Y"))) return(NULL)
  dx <- tr$tip.label[di[m$side[di] == "X"]]
  dy <- tr$tip.label[di[m$side[di] == "Y"]]
  a  <- sub("\\.tre$", "", basename(f))

  res <- list()
  for (s in setdiff(unique(m$genome), c("Nepenthes_gracilis", "Dionaea_muscipula"))) {
    ii <- which(m$genome == s)
    if (length(ii) < 2) next
    for (p in combn(ii, 2, simplify = FALSE)) {
      if (m$chr[p[1]] == m$chr[p[2]]) next
      t1 <- tr$tip.label[p[1]]; t2 <- tr$tip.label[p[2]]
      st <- tryCatch(keep.tip(tr, c(nep, dx, dy, t1, t2)), error = function(e) NULL)
      if (is.null(st) || Ntip(st) != 5) next
      st <- tryCatch(root(st, nep, resolve.root = TRUE), error = function(e) NULL)
      if (is.null(st)) next
      ing <- setdiff(st$tip.label, nep)
      # per-tree null: fraction of the 6 possible ingroup pairs that are clades
      nullp <- mean(apply(combn(ing, 2), 2, function(q) is.monophyletic(st, q)))
      mono <- is.monophyletic(st, c(dx, dy))
      w1 <- NA_character_
      if (!mono) {
        nd <- node.depth(st)
        d1x <- nd[getMRCA(st, c(t1, dx))]; d1y <- nd[getMRCA(st, c(t1, dy))]
        w1 <- if (d1x < d1y) "X" else if (d1x > d1y) "Y" else NA_character_
      }
      res[[length(res) + 1]] <- data.frame(
        anchor = a, sp = s, chr1 = m$chr[p[1]], chr2 = m$chr[p[2]],
        mono = mono, nullp = nullp, chr1_with = w1, stringsAsFactors = FALSE)
    }
  }
  bind_rows(res)
}

q <- bind_rows(lapply(fs, one))
if (!nrow(q)) stop("no quartets extracted")
q <- left_join(q, hits, by = c("anchor" = "nep_gene"))
write_csv(q, "quartet_topology.csv")
cat(sprintf("quartets: %d from %d trees\n", nrow(q), n_distinct(q$anchor)))

cat("\n=== Dionaea monophyly vs per-tree permutation null ===\n")
s <- q %>% group_by(sp) %>%
  summarise(n = n(), obs = sum(mono), exp = sum(nullp),
            frac_obs = mean(mono), frac_null = mean(nullp),
            z = (sum(mono) - sum(nullp)) / sqrt(sum(nullp * (1 - nullp))), .groups = "drop") %>%
  mutate(p = 2 * pnorm(-abs(z)))
print(as.data.frame(s), digits = 3)
cat("frac_obs >> frac_null => Dionaea WGD independent of that species\n")
cat("frac_obs <= frac_null => interleaved => SHARED WGD\n")

cat("\n=== interleaved quartets, by ancestral region ===\n")
print(as.data.frame(q %>% filter(!mono, !is.na(chr1_with)) %>%
  count(sp, nep_chr) %>% pivot_wider(names_from = sp, values_from = n, values_fill = 0)))

## chromosome -> subgenome assignment from interleaved quartets only
assign <- bind_rows(
  q %>% filter(!mono, !is.na(chr1_with)) %>% transmute(sp, chr = chr1, side = chr1_with),
  q %>% filter(!mono, !is.na(chr1_with)) %>%
    transmute(sp, chr = chr2, side = ifelse(chr1_with == "X", "Y", "X"))) %>%
  count(sp, chr, side) %>% pivot_wider(names_from = side, values_from = n, values_fill = 0) %>%
  mutate(tot = X + Y, frac_X = X / tot) %>% filter(tot >= 20)

glob <- assign %>% group_by(sp) %>%
  summarise(base = sum(X) / sum(tot), .groups = "drop")
assign <- assign %>% left_join(glob, by = "sp") %>% rowwise() %>%
  mutate(p_vs_base = binom.test(X, tot, base)$p.value) %>% ungroup() %>%
  group_by(sp) %>% mutate(p_adj = p.adjust(p_vs_base, "BH")) %>% ungroup() %>%
  arrange(sp, desc(frac_X))

cat("\n=== chromosome -> Dionaea subgenome (tested vs each species' own baseline) ===\n")
for (sp in unique(assign$sp)) {
  d <- filter(assign, sp == !!sp)
  cat(sprintf("\n-- %s (baseline frac_X = %.3f) --\n", sp, d$base[1]))
  print(as.data.frame(select(d, chr, X, Y, tot, frac_X, p_adj)), digits = 3)
}
write_csv(assign, "chromosome_subgenome_assignment.csv")
