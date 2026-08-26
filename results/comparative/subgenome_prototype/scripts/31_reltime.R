#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); BAND <- 0.08
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"
SPP <- c("Drosera_binata","Drosera_paradoxa","Drosera_scorpioides",
         "Drosera_regia","Drosera_capensis")

k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS>0, dS<5, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
key <- tm %>% select(anchor=nep_gene, k1=tip, genome, chr) %>% distinct() %>%
       mutate(region=sub("-.*$","",anchor))
DD <- bind_rows(k %>% transmute(anchor, ga=seq1, gb=seq2, d=dS),
                k %>% transmute(anchor, ga=seq2, gb=seq1, d=dS))
nep <- k %>% filter(xor(sp1==NEP, sp2==NEP)) %>%
       transmute(anchor, N=ifelse(sp1==NEP, seq1, seq2)) %>% distinct()
seg <- read_csv("trees/qc/segment_delta_calibrated.csv", show_col_types=FALSE) %>%
  mutate(side=case_when(med_delta < -BAND ~ "A", med_delta > BAND ~ "B", TRUE ~ "amb")) %>%
  select(genome, region, chr, side)

hr("0. build within-species copy pairs, classified by subgenome side")
P <- k %>% filter(sp1==sp2, sp1 %in% SPP) %>%
  transmute(anchor, sp=sp1, g1=seq1, g2=seq2, dSS=dS) %>%
  left_join(key %>% select(anchor, g1=k1, c1=chr, region), by=c("anchor","g1")) %>%
  left_join(key %>% select(anchor, g2=k1, c2=chr), by=c("anchor","g2")) %>%
  filter(!is.na(c1), !is.na(c2), c1 != c2) %>%
  left_join(seg %>% rename(sp=genome, c1=chr, s1=side), by=c("sp","region","c1")) %>%
  left_join(seg %>% rename(sp=genome, c2=chr, s2=side), by=c("sp","region","c2")) %>%
  filter(!is.na(s1), !is.na(s2), s1!="amb", s2!="amb") %>%
  mutate(cls = ifelse(s1==s2, paste0("same-",s1), "cross-AB"))
print(as.data.frame(P %>% count(sp, cls) %>%
  pivot_wider(names_from=cls, values_from=n, values_fill=0)), row.names=FALSE)

hr("1. four-point: are a species' two copies SISTERS relative to another species?")
cat("  sisters  -> the duplication POSTDATES the split from that species\n")
cat("  interleaved -> it PREDATES it\n")
cat("  cross-AB pairs are the control: they must always be interleaved.\n\n")
tg <- bind_rows(k %>% transmute(anchor, g=seq1, tsp=sp1),
                k %>% transmute(anchor, g=seq2, tsp=sp2)) %>%
      distinct() %>% filter(tsp %in% SPP)
Q <- P %>% inner_join(nep, by="anchor") %>%
  inner_join(tg %>% rename(Tg=g), by="anchor", relationship="many-to-many") %>%
  filter(tsp != sp) %>%
  left_join(DD %>% rename(g1=ga, Tg=gb, d_1T=d), by=c("anchor","g1","Tg")) %>%
  left_join(DD %>% rename(g2=ga, Tg=gb, d_2T=d), by=c("anchor","g2","Tg")) %>%
  left_join(DD %>% rename(g1=ga, N=gb,  d_1N=d), by=c("anchor","g1","N")) %>%
  left_join(DD %>% rename(g2=ga, N=gb,  d_2N=d), by=c("anchor","g2","N")) %>%
  left_join(DD %>% rename(Tg=ga, N=gb,  d_TN=d), by=c("anchor","Tg","N")) %>%
  filter(!is.na(d_1T), !is.na(d_2T), !is.na(d_1N), !is.na(d_2N), !is.na(d_TN)) %>%
  mutate(sA = dSS + d_TN, sB = d_1T + d_2N, sC = d_2T + d_1N,
         sisters = (sA <= sB & sA <= sC))
print(as.data.frame(Q %>% group_by(sp, cls) %>%
  summarise(quartets=n(), loci=n_distinct(anchor),
            frac_sisters=round(mean(sisters),3), .groups="drop") %>%
  filter(quartets >= 40) %>% arrange(sp, cls)), row.names=FALSE)
cat("\n  same-A well above the cross-AB control = that A1/A2 split is younger\n")
cat("  than the speciation events, i.e. lineage-specific.\n")

hr("2. same-A only, broken out by which species is the outgroup")
print(as.data.frame(Q %>% filter(grepl("^same-", cls)) %>%
  group_by(sp, tsp) %>%
  summarise(quartets=n(), frac_sisters=round(mean(sisters),3), .groups="drop") %>%
  filter(quartets >= 30) %>%
  pivot_wider(names_from=tsp, values_from=frac_sisters,
              names_prefix="vs_", id_cols=sp)), row.names=FALSE)
cat("\n  High against every outgroup = duplication after all those splits.\n")

hr("3. the control, stated explicitly")
ctl <- Q %>% filter(cls=="cross-AB") %>% group_by(sp) %>%
  summarise(quartets=n(), frac_sisters=round(mean(sisters),3), .groups="drop")
print(as.data.frame(ctl), row.names=FALSE)
cat("\n  A/B homeologs should almost never be sisters (the A/B split is the\n")
cat("  deepest event). A high value here means the side labels are wrong.\n")

hr("4. within-species ratio table — the rate-free depth comparison")
print(as.data.frame(P %>% group_by(sp, cls) %>%
  summarise(pairs=n(), med=round(median(dSS),3), .groups="drop") %>%
  pivot_wider(names_from=cls, values_from=c(pairs, med)) %>%
  mutate(ratio_A = round(`med_same-A`/`med_cross-AB`, 3))), row.names=FALSE)
cat("\n  Both terms are within-species so lineage rate cancels; the two are of\n")
cat("  similar magnitude so saturation largely cancels too.\n")
write_csv(Q %>% select(anchor, sp, cls, tsp, dSS, sA, sB, sC, sisters),
          "trees/qc/reltime_quartets.csv")
