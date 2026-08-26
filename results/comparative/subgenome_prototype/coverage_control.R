#!/usr/bin/env Rscript
# PIN LIBPATH: genespace module sets R_LIBS_USER (rlang 1.1.3, ggplot2 3.5.1)
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

blk <- read_csv(file.path(GSD, "syntenicBlock_coordinates.csv"), show_col_types = FALSE)
bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
frac <- read_csv("fractionation_by_chrpair.csv", show_col_types = FALSE)

## normalise blocks to (nep_chr, nep_start, nep_end, dio_chr)
nd <- blk %>%
  filter((genome1 == "Nepenthes_gracilis" & genome2 == "Dionaea_muscipula") |
         (genome1 == "Dionaea_muscipula"  & genome2 == "Nepenthes_gracilis")) %>%
  mutate(n1 = genome1 == "Nepenthes_gracilis",
         nep_chr = ifelse(n1, chr1, chr2),
         nep_s   = ifelse(n1, pmin(startBp1, endBp1), pmin(startBp2, endBp2)),
         nep_e   = ifelse(n1, pmax(startBp1, endBp1), pmax(startBp2, endBp2)),
         dio_chr = ifelse(n1, chr2, chr1),
         nhit    = ifelse(n1, nHits1, nHits2)) %>%
  filter(grepl("_dom$", nep_chr), !is.na(nep_s), !is.na(nep_e)) %>%
  select(nep_chr, nep_s, nep_e, dio_chr, nhit)

message(sprintf("Nep(dom)-Dio syntenic blocks: %d", nrow(nd)))

merge_iv <- function(s, e) {
  o <- order(s); s <- s[o]; e <- e[o]
  cs <- s[1]; ce <- e[1]; os <- numeric(0); oe <- numeric(0)
  if (length(s) > 1) for (i in 2:length(s)) {
    if (s[i] <= ce) ce <- max(ce, e[i])
    else { os <- c(os, cs); oe <- c(oe, ce); cs <- s[i]; ce <- e[i] }
  }
  list(s = c(os, cs), e = c(oe, ce))
}
in_any <- function(p, s, e) if (!length(s)) rep(FALSE, length(p)) else
  vapply(p, function(x) any(x >= s & x <= e), logical(1))

cov <- nd %>% group_by(nep_chr, dio_chr) %>%
  summarise(n_blk = n(), n_hits = sum(nhit, na.rm = TRUE),
            cov_bp = { m <- merge_iv(nep_s, nep_e); sum(m$e - m$s) }, .groups = "drop")
write_csv(cov, "syntenic_coverage.csv")

## per-pair coverage asymmetry
pc <- frac %>% select(exp_pair, chrA, chrB) %>%
  left_join(select(cov, dio_chr, nep_chr, covA = cov_bp, hitsA = n_hits),
            by = c("chrA" = "dio_chr")) %>%
  left_join(select(cov, dio_chr, nep_chr2 = nep_chr, covB = cov_bp, hitsB = n_hits),
            by = c("chrB" = "dio_chr")) %>%
  filter(nep_chr == nep_chr2) %>%
  mutate(cov_frac_A = covA / (covA + covB), hit_frac_A = hitsA / (hitsA + hitsB)) %>%
  left_join(select(frac, exp_pair, kA, kB, total, frac_A), by = "exp_pair") %>%
  mutate(p_vs_cov = mapply(function(k, t, e) binom.test(k, t, e)$p.value,
                           kA, total, pmin(pmax(cov_frac_A, 1e-6), 1 - 1e-6)),
         p_adj_cov = p.adjust(p_vs_cov, "BH"))

cat("\n=== coverage asymmetry vs retention asymmetry ===\n")
cat("cov_frac_A = chrA's share of Nepenthes bp covered by syntenic blocks.\n")
cat("If frac_A tracks cov_frac_A, the 'fractionation' is a synteny-detection artifact.\n")
print(as.data.frame(select(pc, exp_pair, nep_chr, cov_frac_A, hit_frac_A,
                           frac_A, p_adj_cov)), digits = 3)
ct <- with(pc, cor.test(cov_frac_A, frac_A, method = "spearman"))
cat(sprintf("\nSpearman(cov_frac_A, frac_A) rho=%.3f p=%.4f  <- high rho = artifact\n",
            ct$estimate, ct$p.value))

## retention restricted to DOUBLY-COVERED Nepenthes intervals
npos <- bed %>% filter(genome == "Nepenthes_gracilis") %>%
  transmute(gene = id, nep_chr = chr, nep_pos = (start + end) / 2, globHOG)
dpos <- bed %>% filter(genome == "Dionaea_muscipula") %>%
  mutate(isRep = as.logical(isArrayRep)) %>% filter(is.na(isRep) | isRep) %>%
  transmute(dio_chr = chr, globHOG)

single <- npos %>% filter(grepl("_dom$", nep_chr)) %>%
  group_by(globHOG) %>% filter(n() == 1) %>% ungroup()
dsum <- dpos %>% group_by(globHOG) %>%
  summarise(n_dio = n(), one_chr = ifelse(n() == 1, dio_chr[1], NA), .groups = "drop")
tab <- inner_join(single, dsum, by = "globHOG") %>% filter(n_dio == 1)

res <- lapply(seq_len(nrow(pc)), function(i) {
  r <- pc[i, ]
  ivA <- filter(nd, nep_chr == r$nep_chr, dio_chr == r$chrA)
  ivB <- filter(nd, nep_chr == r$nep_chr, dio_chr == r$chrB)
  if (!nrow(ivA) || !nrow(ivB)) return(NULL)
  mA <- merge_iv(ivA$nep_s, ivA$nep_e); mB <- merge_iv(ivB$nep_s, ivB$nep_e)
  g <- filter(tab, nep_chr == r$nep_chr, one_chr %in% c(r$chrA, r$chrB))
  if (!nrow(g)) return(NULL)
  keep <- in_any(g$nep_pos, mA$s, mA$e) & in_any(g$nep_pos, mB$s, mB$e)
  g <- g[keep, ]
  kA <- sum(g$one_chr == r$chrA); kB <- sum(g$one_chr == r$chrB); tt <- kA + kB
  tibble(exp_pair = r$exp_pair, chrA = r$chrA, chrB = r$chrB,
         kA_dbl = kA, kB_dbl = kB, total_dbl = tt,
         frac_A_dbl = if (tt) kA / tt else NA_real_,
         p_dbl = if (tt >= 5) binom.test(kA, tt, 0.5)$p.value else NA_real_,
         frac_A_all = r$frac_A)
}) %>% bind_rows()

if (nrow(res)) {
  res <- mutate(res, p_adj_dbl = p.adjust(p_dbl, "BH"),
                retained_more_dbl = ifelse(kA_dbl > kB_dbl, chrA, chrB))
  write_csv(res, "fractionation_doubly_covered.csv")
  cat("\n=== retention inside DOUBLY-COVERED intervals only ===\n")
  print(as.data.frame(select(res, exp_pair, kA_dbl, kB_dbl, total_dbl,
                             frac_A_dbl, frac_A_all, p_adj_dbl,
                             retained_more_dbl)), digits = 3)
  cmp <- left_join(res, select(frac, exp_pair, retained_more), by = "exp_pair")
  cat(sprintf("\ndirection preserved in %d/%d pairs after coverage control\n",
              sum(cmp$retained_more_dbl == cmp$retained_more, na.rm = TRUE),
              sum(!is.na(cmp$retained_more_dbl))))
}
