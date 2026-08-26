#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); NPERM <- 500; BAND <- 0.08
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"; DIO <- "Dionaea_muscipula"
CORE <- c("Drosera_binata","Drosera_paradoxa","Drosera_scorpioides")
X <- c("chr1_sg1_s5","chr2_sg1_s3","chr3_sg1_s7","chr4_sg2_s12",
       "chr5_sg2_s16","chr6_sg2_s8","chr7_sg1_s1","chr8_sg1_s4")
spread <- function(v) round(max(v)/min(v), 3)

k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS>0, dS<3, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
key <- tm %>% select(anchor=nep_gene, k1=tip, genome, chr) %>% distinct() %>%
       mutate(region=sub("-.*$","",anchor))
RT <- read_csv("trees/qc/lineage_rates.csv", show_col_types=FALSE) %>%
      filter(dio_est=="dmin") %>% select(sp=species, rate)
DD <- bind_rows(k %>% transmute(anchor, ga=seq1, gb=seq2, d=dS),
                k %>% transmute(anchor, ga=seq2, gb=seq1, d=dS))
SEG <- read_csv("trees/qc/segment_delta_calibrated.csv", show_col_types=FALSE)

hr("ATTACK 5 — is cross-AB doing more than selecting divergent pairs?")
cat("  cross-AB pairs are BY CONSTRUCTION the more divergent ones (opposite\n")
cat("  sides of the axis). If simply taking the most divergent pairs from each\n")
cat("  species collapses as well, the delta axis contributes nothing.\n\n")
seg <- SEG %>% mutate(side=case_when(med_delta < -BAND ~ "A", med_delta > BAND ~ "B",
                                     TRUE ~ "amb")) %>% select(genome, region, chr, side)
W <- k %>% filter(sp1==sp2, sp1!=NEP) %>% transmute(anchor, sp=sp1, g1=seq1, g2=seq2, dS) %>%
  left_join(key %>% select(anchor, g1=k1, c1=chr, region), by=c("anchor","g1")) %>%
  left_join(key %>% select(anchor, g2=k1, c2=chr), by=c("anchor","g2")) %>%
  filter(!is.na(c1), !is.na(c2), c1!=c2) %>% left_join(RT, by="sp")
lab <- W %>%
  left_join(seg %>% rename(sp=genome, c1=chr, s1=side), by=c("sp","region","c1")) %>%
  left_join(seg %>% rename(sp=genome, c2=chr, s2=side), by=c("sp","region","c2"))
NC <- lab %>% filter(!is.na(s1), !is.na(s2), s1!="amb", s2!="amb", s1!=s2) %>%
      count(sp, name="n_cross")
obs <- lab %>% filter(!is.na(s1), !is.na(s2), s1!="amb", s2!="amb", s1!=s2) %>%
  group_by(sp) %>% summarise(T=median(dS)/(2*rate[1]), .groups="drop")
print(as.data.frame(W %>% count(sp, name="n_all") %>% left_join(NC, by="sp") %>%
  mutate(frac=round(n_cross/n_all,3)) %>% left_join(obs %>% rename(T_cross=T), by="sp") %>%
  mutate(T_cross=round(T_cross,3))), row.names=FALSE)
cat(sprintf("\n  observed cross-AB spread: %.3fx\n\n", spread(obs$T)))

topn <- W %>% inner_join(NC, by="sp") %>% group_by(sp) %>%
  slice_max(dS, n=n_cross[1]) %>%
  summarise(T=median(dS)/(2*rate[1]), .groups="drop")
cat(sprintf("  top-N by dS (same N per species) : %s   spread %.3fx\n",
    paste(round(sort(topn$T),3), collapse=" "), spread(topn$T)))
set.seed(1)
rnd <- vapply(seq_len(NPERM), function(i) {
  q <- W %>% inner_join(NC, by="sp") %>% group_by(sp) %>%
       slice_sample(n=n_cross[1]) %>%
       summarise(T=median(dS)/(2*rate[1]), .groups="drop")
  spread(q$T)
}, numeric(1))
cat(sprintf("  random N per species             : median %.3fx  [%.3f - %.3f]\n",
            median(rnd), quantile(rnd,.05), quantile(rnd,.95)))
cat(sprintf("\n  random splits at least as tight as observed: %.3f\n", mean(rnd <= spread(obs$T))))
cat("  Top-N also collapsing = the axis is redundant with divergence.\n")
cat("  Top-N staying spread = the axis is selecting something else.\n")

hr("ATTACK 6 — is the 2:1 constitution robust to a random orientation?")
cat("  This is the version of attack 2 that is NOT a no-op. Flipping a region\n")
cat("  swaps its A and B counts, so under the chromosome-autonomous null\n")
cat("  (v2 sec0 case b) each region reads 2A:1B or 1A:2B at random and the\n")
cat("  totals should average to 1:1.\n\n")
obs_c <- seg %>% filter(side!="amb", genome %in% CORE) %>% group_by(genome) %>%
  summarise(A=sum(side=="A"), B=sum(side=="B"),
            ratio=round(sum(side=="A")/pmax(1,sum(side=="B")),3), .groups="drop")
print(as.data.frame(obs_c), row.names=FALSE)
regions <- sort(unique(seg$region))
set.seed(2)
nulc <- lapply(seq_len(NPERM), function(i) {
  fl <- sample(regions, rbinom(1, length(regions), 0.5))
  seg %>% filter(side!="amb", genome %in% CORE) %>%
    mutate(side=ifelse(region %in% fl, ifelse(side=="A","B","A"), side)) %>%
    group_by(genome) %>%
    summarise(ratio=sum(side=="A")/pmax(1,sum(side=="B")), .groups="drop")
})
nd <- bind_rows(nulc)
cat("\n  under random per-region orientation:\n")
print(as.data.frame(nd %>% group_by(genome) %>%
  summarise(med=round(median(ratio),3), lo=round(quantile(ratio,.025),3),
            hi=round(quantile(ratio,.975),3), .groups="drop") %>%
  left_join(obs_c %>% select(genome, observed=ratio), by="genome") %>%
  rowwise() %>%
  mutate(p=round(mean(nd$ratio[nd$genome==genome] >= observed),4)) %>% ungroup()),
  row.names=FALSE)
cat("\n  Observed far above the null = the fractionation orientation IS globally\n")
cat("  consistent, which is evidence AGAINST v2's chromosome-autonomous case (b)\n")
cat("  and would partly rescue the assumption flagged in the assessment.\n")
cat("  Observed inside the null = the 2:1 could be an arbitrary partition.\n")

hr("ATTACK 7 — is the A-heavy skew itself orientation-dependent?")
cat("  Global share of segments called A, vs the same under random flips.\n\n")
oa <- seg %>% filter(side!="amb") %>% summarise(fA=mean(side=="A")) %>% pull(fA)
na_ <- vapply(seq_len(NPERM), function(i) {
  fl <- sample(regions, rbinom(1, length(regions), 0.5))
  seg %>% filter(side!="amb") %>%
    mutate(side=ifelse(region %in% fl, ifelse(side=="A","B","A"), side)) %>%
    summarise(f=mean(side=="A")) %>% pull(f)
}, numeric(1))
cat(sprintf("  observed A share : %.3f\n", oa))
cat(sprintf("  random flips     : median %.3f  [%.3f - %.3f]\n",
            median(na_), quantile(na_,.025), quantile(na_,.975)))
cat(sprintf("  p (random >= observed): %.4f\n", mean(na_ >= oa)))
