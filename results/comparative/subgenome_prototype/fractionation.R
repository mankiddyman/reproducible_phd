#!/usr/bin/env Rscript
# PIN LIBPATH: genespace module sets R_LIBS_USER (rlang 1.1.3, ggplot2 3.5.1)
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))

# Retention asymmetry = second sign vector over the same 8 pairs, mechanistically
# independent of dS. Nep->Dio pair mapping calibrated from 1:2 HOGs, not assumed.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

GS <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results/combBed.txt"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

bed <- read.table(GS, header = TRUE, sep = "\t", quote = "", comment.char = "",
                  stringsAsFactors = FALSE)

b <- bed %>%
  filter(genome %in% c("Nepenthes_gracilis", "Dionaea_muscipula")) %>%
  mutate(isRep = as.logical(isArrayRep)) %>%
  filter(is.na(isRep) | isRep) %>%                 # collapse tandem arrays
  mutate(pair_lab = ifelse(genome == "Dionaea_muscipula",
                           sub("_sg[0-9]+_s[0-9]+$", "", chr), NA_character_))

message(sprintf("array reps: Nep %d | Dio %d",
                sum(b$genome == "Nepenthes_gracilis"), sum(b$genome == "Dionaea_muscipula")))

nep1 <- b %>% filter(genome == "Nepenthes_gracilis", grepl("_dom$", chr)) %>%
  group_by(globHOG) %>% filter(n() == 1) %>% ungroup() %>%
  select(globHOG, nep_chr = chr)

dio <- b %>% filter(genome == "Dionaea_muscipula") %>%
  group_by(globHOG) %>%
  summarise(n_dio = n(),
            pairs = paste(sort(unique(pair_lab)), collapse = ";"),
            chrs  = paste(sort(chr), collapse = ";"),
            one_chr = ifelse(n() == 1, chr[1], NA_character_), .groups = "drop")

tab <- inner_join(nep1, dio, by = "globHOG")

# calibrate: which Dionaea pair does each Nepenthes dom chromosome map to?
cal <- tab %>% filter(n_dio == 2, !grepl(";", pairs)) %>%
  count(nep_chr, pairs, name = "cal_n") %>%
  group_by(nep_chr) %>% mutate(frac = cal_n / sum(cal_n)) %>%
  slice_max(cal_n, n = 1, with_ties = FALSE) %>% ungroup() %>%
  rename(exp_pair = pairs)
cat("\n=== Nep chr -> Dio pair calibration (frac = purity) ===\n")
print(as.data.frame(cal), digits = 3)
stopifnot(!any(duplicated(cal$exp_pair)))   # mapping must be 1:1

tab  <- left_join(tab, select(cal, nep_chr, exp_pair), by = "nep_chr")
lost <- tab %>% filter(n_dio == 1, !is.na(exp_pair),
                       sub("_sg[0-9]+_s[0-9]+$", "", one_chr) == exp_pair)
retn <- tab %>% filter(n_dio == 2, !is.na(exp_pair), pairs == exp_pair)

frac <- lost %>% count(exp_pair, one_chr, name = "k") %>%
  group_by(exp_pair) %>%
  complete(one_chr, fill = list(k = 0L)) %>%
  arrange(exp_pair, one_chr) %>%
  summarise(chrA = one_chr[1], chrB = one_chr[2], kA = k[1], kB = k[2], .groups = "drop") %>%
  mutate(total = kA + kB, frac_A = kA / total,
         p = mapply(function(a, t) if (t >= 5) binom.test(a, t, 0.5)$p.value else NA_real_,
                    kA, total),
         p_adj = p.adjust(p, "BH"),
         retained_more = ifelse(kA > kB, chrA, chrB)) %>%
  left_join(count(retn, exp_pair, name = "n_1to2"), by = "exp_pair")

write_csv(frac, "fractionation_by_chrpair.csv")

if (file.exists("phasing_by_chrpair_v2.csv")) {
  ds <- read_csv("phasing_by_chrpair_v2.csv", show_col_types = FALSE) %>%
    filter(set == "conv_filtered") %>% select(pair_lab, ds_faster = faster,
                                              ds_est = est, ds_n = n)
  cc <- frac %>% left_join(ds, by = c("exp_pair" = "pair_lab")) %>%
    # dominance predicts faster-evolving == more-fractionated == kept FEWER
    mutate(frac_loser = ifelse(kA > kB, chrB, chrA),
           agree = ds_faster == frac_loser)
  write_csv(cc, "concordance_ds_vs_fractionation.csv")
  k <- sum(cc$agree, na.rm = TRUE); nn <- sum(!is.na(cc$agree))
  cat("\n=== concordance: dS-faster vs more-fractionated ===\n")
  print(as.data.frame(select(cc, exp_pair, chrA, chrB, kA, kB, retained_more,
                             ds_faster, ds_est, agree)), digits = 3)
  if (nn > 0) cat(sprintf("\nagree %d/%d | binomial p=%.4f | floor with one label fixed = %.4f\n",
                          k, nn, binom.test(max(k, nn - k), nn, 0.5)$p.value, 0.5^(nn - 1)))
}

cat("\n=== fractionation ===\n"); print(as.data.frame(frac), digits = 3)
cat(sprintf("\nHOGs Nep-single-on-dom %d | 1:1 usable %d | 1:2 usable %d\n",
            nrow(tab), nrow(lost), nrow(retn)))
