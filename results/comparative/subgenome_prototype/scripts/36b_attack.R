#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); NPERM <- 500; BAND <- 0.08
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"; DIO <- "Dionaea_muscipula"
CORE <- c("Drosera_binata","Drosera_paradoxa","Drosera_scorpioides")
spread <- function(v) round(max(v)/min(v), 3)

k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS>0, dS<3, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
key <- tm %>% select(anchor=nep_gene, k1=tip, genome, chr) %>% distinct() %>%
       mutate(region=sub("-.*$","",anchor))
RT <- read_csv("trees/qc/lineage_rates.csv", show_col_types=FALSE) %>%
      filter(dio_est=="dmin") %>% select(sp=species, rate)
SEG <- read_csv("trees/qc/segment_delta_calibrated.csv", show_col_types=FALSE)
seg <- SEG %>% mutate(side=case_when(med_delta < -BAND ~ "A", med_delta > BAND ~ "B",
                                     TRUE ~ "amb")) %>% select(genome, region, chr, side)
W <- k %>% filter(sp1==sp2, sp1!=NEP, sp1!=DIO) %>%
  transmute(anchor, sp=sp1, g1=seq1, g2=seq2, dS) %>%
  left_join(key %>% select(anchor, g1=k1, c1=chr, region), by=c("anchor","g1")) %>%
  left_join(key %>% select(anchor, g2=k1, c2=chr), by=c("anchor","g2")) %>%
  filter(!is.na(c1), !is.na(c2), c1!=c2) %>% left_join(RT, by="sp")
lab <- W %>%
  left_join(seg %>% rename(sp=genome, c1=chr, s1=side), by=c("sp","region","c1")) %>%
  left_join(seg %>% rename(sp=genome, c2=chr, s2=side), by=c("sp","region","c2"))
CR <- lab %>% filter(!is.na(s1), !is.na(s2), s1!="amb", s2!="amb", s1!=s2)
NC <- CR %>% count(sp, name="N")
obs <- CR %>% group_by(sp) %>% summarise(T=median(dS)/(2*rate[1]), .groups="drop")
Tall <- W %>% group_by(sp) %>% summarise(T=median(dS)/(2*rate[1]), .groups="drop")

hr("ATTACK 5 (fixed) — does a size-matched subset collapse as well?")
pick <- function(d, mode) { N <- d$N[1]
  if (mode=="top") d[order(-d$dS),][seq_len(N),]
  else if (mode=="bot") d[order(d$dS),][seq_len(N),]
  else d[sample(nrow(d), N),] }
sub_T <- function(mode) W %>% inner_join(NC, by="sp") %>% group_by(sp) %>%
  group_modify(~pick(.x, mode)) %>%
  summarise(T=median(dS)/(2*rate[1]), .groups="drop")
cat(sprintf("  all pairs, no labels : %s   spread %.3fx\n",
    paste(round(sort(Tall$T),3), collapse=" "), spread(Tall$T)))
cat(sprintf("  cross-AB (observed)  : %s   spread %.3fx\n",
    paste(round(sort(obs$T),3), collapse=" "), spread(obs$T)))
tp <- sub_T("top"); bt <- sub_T("bot")
cat(sprintf("  top-N by dS          : %s   spread %.3fx\n",
    paste(round(sort(tp$T),3), collapse=" "), spread(tp$T)))
cat(sprintf("  bottom-N by dS       : %s   spread %.3fx\n",
    paste(round(sort(bt$T),3), collapse=" "), spread(bt$T)))
set.seed(1)
rnd <- vapply(seq_len(NPERM), function(i) spread(sub_T("rnd")$T), numeric(1))
cat(sprintf("  random N per species : median %.3fx  [%.3f - %.3f]\n",
            median(rnd), quantile(rnd,.05), quantile(rnd,.95)))
cat(sprintf("\n  random subsets at least as tight as observed: %.4f\n",
            mean(rnd <= spread(obs$T))))
cat("\n  excluding capensis:\n")
nc <- function(d) d %>% filter(sp!="Drosera_capensis")
cat(sprintf("    all pairs %.3fx | cross-AB %.3fx | random median %.3fx\n",
  spread(nc(Tall)$T), spread(nc(obs)$T),
  median(vapply(seq_len(NPERM), function(i) spread(nc(sub_T("rnd"))$T), numeric(1)))))

hr("ATTACK 6 — is the 2:1 constitution robust to a random orientation?")
cat("  Under v2 sec0 case (b), X/Y is arbitrary w.r.t. ancestry, so each\n")
cat("  region reads 2A:1B or 1A:2B at random and totals average to 1:1.\n\n")
obs_c <- seg %>% filter(side!="amb", genome %in% CORE) %>% group_by(genome) %>%
  summarise(A=sum(side=="A"), B=sum(side=="B"),
            ratio=round(sum(side=="A")/pmax(1,sum(side=="B")),3), .groups="drop")
print(as.data.frame(obs_c), row.names=FALSE)
regions <- sort(unique(seg$region))
set.seed(2)
nd <- bind_rows(lapply(seq_len(NPERM), function(i) {
  fl <- sample(regions, rbinom(1, length(regions), 0.5))
  seg %>% filter(side!="amb", genome %in% CORE) %>%
    mutate(side=ifelse(region %in% fl, ifelse(side=="A","B","A"), side)) %>%
    group_by(genome) %>%
    summarise(ratio=sum(side=="A")/pmax(1,sum(side=="B")), .groups="drop") }))
cat("\n  under random per-region orientation:\n")
print(as.data.frame(nd %>% group_by(genome) %>%
  summarise(med=round(median(ratio),3), lo=round(quantile(ratio,.025),3),
            hi=round(quantile(ratio,.975),3), .groups="drop") %>%
  left_join(obs_c %>% select(genome, observed=ratio), by="genome") %>%
  rowwise() %>% mutate(p=round(mean(nd$ratio[nd$genome==genome] >= observed),4)) %>%
  ungroup()), row.names=FALSE)
cat("\n  Observed above the null = the orientation IS globally consistent,\n")
cat("  which is evidence against case (b). Inside the null = arbitrary.\n")

hr("ATTACK 7 — is the global A-heavy skew orientation-dependent?")
oa <- seg %>% filter(side!="amb") %>% summarise(f=mean(side=="A")) %>% pull(f)
na_ <- vapply(seq_len(NPERM), function(i) {
  fl <- sample(regions, rbinom(1, length(regions), 0.5))
  seg %>% filter(side!="amb") %>%
    mutate(side=ifelse(region %in% fl, ifelse(side=="A","B","A"), side)) %>%
    summarise(f=mean(side=="A")) %>% pull(f) }, numeric(1))
cat(sprintf("  observed A share %.3f | random median %.3f [%.3f - %.3f] | p = %.4f\n",
            oa, median(na_), quantile(na_,.025), quantile(na_,.975), mean(na_ >= oa)))
