#!/usr/bin/env Rscript
# PIN LIBPATH: genespace module sets R_LIBS_USER (rlang 1.1.3, ggplot2 3.5.1)
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Same fractionation test, two independent Nepenthes anchor sets (_dom vs _nondom).
# Agreement => partition is a Dionaea property, not an anchor artifact.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
blk <- read_csv(file.path(GSD, "syntenicBlock_coordinates.csv"), show_col_types = FALSE)

cat("=== Nepenthes chromosome inventory ===\n")
print(bed %>% filter(genome == "Nepenthes_gracilis") %>% count(chr) %>% as.data.frame())

b <- bed %>% filter(genome %in% c("Nepenthes_gracilis", "Dionaea_muscipula")) %>%
  mutate(isRep = as.logical(isArrayRep)) %>% filter(is.na(isRep) | isRep)

nd_all <- blk %>%
  filter((genome1 == "Nepenthes_gracilis" & genome2 == "Dionaea_muscipula") |
         (genome1 == "Dionaea_muscipula"  & genome2 == "Nepenthes_gracilis")) %>%
  mutate(n1 = genome1 == "Nepenthes_gracilis",
         nep_chr = ifelse(n1, chr1, chr2),
         nep_s = ifelse(n1, pmin(startBp1, endBp1), pmin(startBp2, endBp2)),
         nep_e = ifelse(n1, pmax(startBp1, endBp1), pmax(startBp2, endBp2)),
         dio_chr = ifelse(n1, chr2, chr1)) %>%
  filter(!is.na(nep_s), !is.na(nep_e)) %>% select(nep_chr, nep_s, nep_e, dio_chr)

merge_iv <- function(s, e) {
  o <- order(s); s <- s[o]; e <- e[o]; cs <- s[1]; ce <- e[1]
  os <- numeric(0); oe <- numeric(0)
  if (length(s) > 1) for (i in 2:length(s)) {
    if (s[i] <= ce) ce <- max(ce, e[i]) else { os <- c(os, cs); oe <- c(oe, ce); cs <- s[i]; ce <- e[i] }
  }
  list(s = c(os, cs), e = c(oe, ce))
}
in_any <- function(p, s, e) if (!length(s)) rep(FALSE, length(p)) else
  vapply(p, function(x) any(x >= s & x <= e), logical(1))

run_anchor <- function(rx, tag) {
  nep <- b %>% filter(genome == "Nepenthes_gracilis", grepl(rx, chr)) %>%
    transmute(globHOG, nep_chr = chr, nep_pos = (start + end)/2)
  if (!nrow(nep)) { message("no Nepenthes chromosomes match ", rx); return(NULL) }
  nep1 <- nep %>% group_by(globHOG) %>% filter(n() == 1) %>% ungroup()

  dio <- b %>% filter(genome == "Dionaea_muscipula") %>%
    mutate(pl = sub("_sg[0-9]+_s[0-9]+$", "", chr)) %>%
    group_by(globHOG) %>%
    summarise(n_dio = n(), pairs = paste(sort(unique(pl)), collapse = ";"),
              one_chr = ifelse(n() == 1, chr[1], NA), .groups = "drop")

  tab <- inner_join(nep1, dio, by = "globHOG")
  cal <- tab %>% filter(n_dio == 2, !grepl(";", pairs)) %>%
    count(nep_chr, pairs, name = "cal_n") %>% group_by(nep_chr) %>%
    mutate(purity = cal_n/sum(cal_n)) %>% slice_max(cal_n, n = 1, with_ties = FALSE) %>%
    ungroup() %>% rename(exp_pair = pairs)
  cat(sprintf("\n=== [%s] calibration ===\n", tag)); print(as.data.frame(cal), digits = 3)
  if (any(duplicated(cal$exp_pair))) cat("!! non-bijective mapping — interpret with care\n")

  nd <- filter(nd_all, grepl(rx, nep_chr))
  lost <- tab %>% inner_join(select(cal, nep_chr, exp_pair), by = "nep_chr") %>%
    filter(n_dio == 1, sub("_sg[0-9]+_s[0-9]+$", "", one_chr) == exp_pair)

  out <- lapply(split(lost, lost$exp_pair), function(g) {
    chrs <- sort(unique(sub("^.*$", "", character(0))))
    ch <- sort(unique(g$one_chr))
    if (length(ch) != 2) return(NULL)
    nc <- g$nep_chr[1]
    ivA <- filter(nd, nep_chr == nc, dio_chr == ch[1])
    ivB <- filter(nd, nep_chr == nc, dio_chr == ch[2])
    kA <- sum(g$one_chr == ch[1]); kB <- sum(g$one_chr == ch[2])
    dblA <- dblB <- NA_integer_
    if (nrow(ivA) && nrow(ivB)) {
      mA <- merge_iv(ivA$nep_s, ivA$nep_e); mB <- merge_iv(ivB$nep_s, ivB$nep_e)
      k <- in_any(g$nep_pos, mA$s, mA$e) & in_any(g$nep_pos, mB$s, mB$e)
      dblA <- sum(g$one_chr[k] == ch[1]); dblB <- sum(g$one_chr[k] == ch[2])
    }
    tibble(anchor = tag, exp_pair = g$exp_pair[1], nep_chr = nc,
           chrA = ch[1], chrB = ch[2], kA = kA, kB = kB,
           frac_A = kA/(kA+kB), kA_dbl = dblA, kB_dbl = dblB,
           frac_A_dbl = if (is.na(dblA) || (dblA+dblB) == 0) NA_real_ else dblA/(dblA+dblB),
           winner = ifelse(kA > kB, ch[1], ch[2]),
           winner_dbl = if (is.na(dblA)) NA_character_ else ifelse(dblA > dblB, ch[1], ch[2]))
  }) %>% bind_rows()

  if (nrow(out)) out <- mutate(out,
    p = mapply(function(a,bb) if (a+bb >= 5) binom.test(a, a+bb, 0.5)$p.value else NA_real_, kA, kB),
    p_adj = p.adjust(p, "BH"))
  cat(sprintf("\n=== [%s] fractionation ===\n", tag))
  print(as.data.frame(select(out, exp_pair, nep_chr, kA, kB, frac_A,
                             frac_A_dbl, p_adj, winner, winner_dbl)), digits = 3)
  out
}

dom <- run_anchor("_dom$", "dom")
non <- run_anchor("_nondom$", "nondom")

if (!is.null(dom) && !is.null(non)) {
  cmp <- inner_join(select(dom, exp_pair, w_dom = winner, wd_dom = winner_dbl,
                           f_dom = frac_A, chrA, chrB),
                    select(non, exp_pair, w_non = winner, wd_non = winner_dbl,
                           f_non = frac_A), by = "exp_pair") %>%
    mutate(agree = w_dom == w_non, agree_dbl = wd_dom == wd_non)
  write_csv(cmp, "anchor_replication.csv")
  cat("\n=== dom vs nondom anchor ===\n"); print(as.data.frame(cmp), digits = 3)
  k <- sum(cmp$agree, na.rm = TRUE); n <- sum(!is.na(cmp$agree))
  cat(sprintf("\nagree %d/%d | sign test p=%.4f\n", k, n,
              binom.test(max(k, n-k), n, 0.5)$p.value))
  kd <- sum(cmp$agree_dbl, na.rm = TRUE); nd2 <- sum(!is.na(cmp$agree_dbl))
  if (nd2) cat(sprintf("doubly-covered: agree %d/%d | p=%.4f\n", kd, nd2,
                       binom.test(max(kd, nd2-kd), nd2, 0.5)$p.value))
  ct <- cor.test(cmp$f_dom, cmp$f_non, method = "spearman")
  cat(sprintf("Spearman(frac_dom, frac_nondom) rho=%.3f p=%.4f\n", ct$estimate, ct$p.value))
}
