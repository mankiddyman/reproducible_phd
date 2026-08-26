#!/usr/bin/env Rscript
# PIN LIBPATH: genespace module sets R_LIBS_USER (rlang 1.1.3, ggplot2 3.5.1)
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
FAI <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/Dionaea_muscipula/assembly_final/external_collapsed/Dionaea_muscipula_chr.fa.fai"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

bed  <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                   quote = "", comment.char = "", stringsAsFactors = FALSE)
frac <- read_csv("fractionation_by_chrpair.csv", show_col_types = FALSE)
ds   <- read_csv("phasing_by_chrpair_v2.csv", show_col_types = FALSE) %>% filter(set == "conv_filtered")
fai  <- read.table(FAI, stringsAsFactors = FALSE)[,1:2]; names(fai) <- c("chr","len")

## ---- CONTROL 1: total gene content per chromosome ---------------------------
gc <- bed %>% filter(genome == "Dionaea_muscipula") %>%
  mutate(isRep = as.logical(isArrayRep)) %>% filter(is.na(isRep) | isRep) %>%
  count(chr, name = "n_genes") %>% left_join(fai, by = "chr") %>%
  mutate(dens = n_genes / (len/1e6))

ct1 <- frac %>%
  left_join(select(gc, chr, gA = n_genes, lA = len, dA = dens), by = c("chrA"="chr")) %>%
  left_join(select(gc, chr, gB = n_genes, lB = len, dB = dens), by = c("chrB"="chr")) %>%
  mutate(exp_A  = gA/(gA+gB),
         obs_A  = frac_A,
         excess = obs_A - exp_A,
         p_vs_content = mapply(function(k,t,e) binom.test(k,t,e)$p.value, kA, total, exp_A),
         p_adj_content = p.adjust(p_vs_content, "BH"))
cat("\n=== CONTROL 1: retention vs total gene content ===\n")
cat("exp_A = chrA's share of annotated genes. excess>0 = fractionation beyond content asymmetry.\n")
print(as.data.frame(select(ct1, exp_pair, chrA, chrB, gA, gB, exp_A, obs_A,
                           excess, p_adj_content)), digits = 3)

## ---- CONTROL 2: gene density, X vs Y across pairs (cross-pair coherence) ----
part <- frac %>% transmute(exp_pair, X = retained_more,
                           Y = ifelse(retained_more == chrA, chrB, chrA))
coh <- part %>%
  left_join(select(gc, chr, densX = dens, lenX = len), by = c("X"="chr")) %>%
  left_join(select(gc, chr, densY = dens, lenY = len), by = c("Y"="chr")) %>%
  mutate(d_dens = densX - densY, d_len_Mb = (lenX - lenY)/1e6)
cat("\n=== CONTROL 2: does the fractionation partition show coherent gene density? ===\n")
print(as.data.frame(coh), digits = 3)
cat(sprintf("density X vs Y, paired Wilcoxon p=%.3f | sign %d/8 positive\n",
            wilcox.test(coh$densX, coh$densY, paired = TRUE)$p.value, sum(coh$d_dens > 0)))
cat(sprintf("length  X vs Y, paired Wilcoxon p=%.3f | sign %d/8 positive\n",
            wilcox.test(coh$lenX, coh$lenY, paired = TRUE)$p.value, sum(coh$d_len_Mb > 0)))

## ---- CONTROL 3: direction-agnostic concordance ------------------------------
cc <- frac %>% left_join(select(ds, pair_lab, ds_faster = faster, ds_est = est,
                                ds_lo = lo, ds_hi = hi, ds_n = n),
                         by = c("exp_pair" = "pair_lab")) %>%
  mutate(faster_is_retained = ds_faster == retained_more,
         ds_informative = (ds_lo > 0) | (ds_hi < 0))
k <- sum(cc$faster_is_retained); n <- nrow(cc)
ki <- sum(cc$faster_is_retained[cc$ds_informative]); ni <- sum(cc$ds_informative)
cat("\n=== CONTROL 3: concordance, direction-agnostic ===\n")
print(as.data.frame(select(cc, exp_pair, retained_more, frac_A, ds_faster, ds_est,
                           ds_n, ds_informative, faster_is_retained)), digits = 3)
cat(sprintf("\nfaster == more-retained in %d/%d pairs | two-tailed sign test p=%.4f\n",
            k, n, binom.test(max(k, n-k), n, 0.5)$p.value))
cat(sprintf("restricted to pairs whose dS CI excludes zero: %d/%d | p=%.4f\n",
            ki, ni, if (ni > 0) binom.test(max(ki, ni-ki), ni, 0.5)$p.value else NA))

## ---- CONTROL 4: random-effects pooling of dS --------------------------------
s <- filter(ds, is.finite(est), hi > lo)
se <- (s$hi - s$lo)/(2*qnorm(0.975)); w <- 1/se^2
Q <- sum(w*(s$est - sum(w*s$est)/sum(w))^2); dfq <- nrow(s)-1
tau2 <- max(0, (Q - dfq)/(sum(w) - sum(w^2)/sum(w)))
ws <- 1/(se^2 + tau2); mu <- sum(ws*s$est)/sum(ws); semu <- sqrt(1/sum(ws))
cat(sprintf("\n=== CONTROL 4: dS random-effects pooling ===\ntau=%.5f  pooled=%.5f  SE=%.5f  z=%.2f  p=%.3f\n",
            sqrt(tau2), mu, semu, mu/semu, 2*pnorm(-abs(mu/semu))))

## ---- proposed phasing -------------------------------------------------------
cat("\n=== PROVISIONAL PARTITION (fractionation-defined) ===\n")
cat("X (retains more):", paste(sort(part$X), collapse = " "), "\n")
cat("Y (retains less):", paste(sort(part$Y), collapse = " "), "\n")
write_csv(part, "provisional_partition.csv")
write_csv(cc, "concordance_v2.csv")
