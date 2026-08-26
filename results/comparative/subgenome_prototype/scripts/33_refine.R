#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"; DIO <- "Dionaea_muscipula"

k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS>0, dS<3, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
key <- tm %>% select(anchor=nep_gene, k1=tip, genome, chr) %>% distinct() %>%
       mutate(region=sub("-.*$","",anchor))
SEG <- read_csv("trees/qc/segment_delta_calibrated.csv", show_col_types=FALSE)
RT <- read_csv("trees/qc/lineage_rates.csv", show_col_types=FALSE) %>%
      filter(dio_est=="dmin") %>% select(sp=species, rate)
cat("rates in use:\n"); print(as.data.frame(RT), row.names=FALSE)

sideof <- function(B) SEG %>%
  mutate(side=case_when(med_delta < -B ~ "A", med_delta > B ~ "B", TRUE ~ "amb")) %>%
  filter(side != "amb") %>% select(genome, region, chr, side)

hr("1. BAND SWEEP — does same-A converge as contamination is removed?")
cat("  Mislabelled cross-AB pairs inflate the same-A median toward the A/B\n")
cat("  value, so tightening the band should pull same-A DOWN.\n\n")
sweep <- bind_rows(lapply(c(0.08, 0.15, 0.20, 0.25), function(B) {
  sg <- sideof(B)
  k %>% filter(sp1==sp2, sp1!=NEP) %>% transmute(anchor, sp=sp1, g1=seq1, g2=seq2, dS) %>%
    left_join(key %>% select(anchor, g1=k1, c1=chr, region), by=c("anchor","g1")) %>%
    left_join(key %>% select(anchor, g2=k1, c2=chr), by=c("anchor","g2")) %>%
    filter(!is.na(c1), !is.na(c2), c1!=c2) %>%
    left_join(sg %>% rename(sp=genome, c1=chr, s1=side), by=c("sp","region","c1")) %>%
    left_join(sg %>% rename(sp=genome, c2=chr, s2=side), by=c("sp","region","c2")) %>%
    filter(!is.na(s1), !is.na(s2)) %>%
    mutate(cls=ifelse(s1==s2, paste0("same-",s1), "cross-AB"), band=B)
}))
res <- sweep %>% left_join(RT, by="sp") %>% group_by(sp, cls, band) %>%
  summarise(n=n(), T=round(median(dS)/(2*rate[1]),3), .groups="drop")
print(as.data.frame(res %>% filter(cls %in% c("cross-AB","same-A")) %>%
  pivot_wider(names_from=band, values_from=c(n,T)) %>%
  arrange(cls, sp)), row.names=FALSE)
cat("\n  cross-AB should stay flat near 0.57 (it is the shared A/B split).\n")
cat("  same-A dropping and the core three converging = contamination confirmed.\n")

hr("2. same-B counts — and why they are only weak evidence")
print(as.data.frame(sweep %>% filter(cls=="same-B") %>% count(sp, band) %>%
  pivot_wider(names_from=band, values_from=n, values_fill=0)), row.names=FALSE)
cat("\n  Zero for the 2A:1B species, but this is partly forced: 2 A labels and\n")
cat("  1 B label per region make a same-B pair impossible by construction.\n")
cat("  The copy-level version is script 29 section 3 (A=0 depletion).\n")

hr("3. A-LINEAGE vs B-LINEAGE divergence times")
cat("  One copy per side per species, so no min-of-many bias on either side.\n")
cat("  If one hybridisation seeded both lineages, T_A must equal T_B.\n\n")
sg <- sideof(0.15)
lab <- key %>% left_join(sg, by=c("genome","region","chr")) %>%
       filter(side %in% c("A","B"))
one <- lab %>% count(anchor, genome, side, name="k") %>% filter(k==1)
CS <- k %>% filter(sp1!=sp2, sp1!=NEP, sp2!=NEP, sp1!=DIO, sp2!=DIO) %>%
  left_join(lab %>% select(anchor, seq1=k1, s1=side), by=c("anchor","seq1")) %>%
  left_join(lab %>% select(anchor, seq2=k1, s2=side), by=c("anchor","seq2")) %>%
  filter(!is.na(s1), !is.na(s2), s1==s2) %>%
  semi_join(one, by=c("anchor","sp1"="genome","s1"="side")) %>%
  semi_join(one, by=c("anchor","sp2"="genome","s2"="side")) %>%
  mutate(a=pmin(sp1,sp2), b=pmax(sp1,sp2), side=s1)
cmp <- CS %>% left_join(RT %>% rename(a=sp, ra=rate), by="a") %>%
  left_join(RT %>% rename(b=sp, rb=rate), by="b") %>%
  group_by(a, b, side) %>%
  summarise(n=n(), T=round(median(dS)/(ra[1]+rb[1]),3), .groups="drop") %>%
  pivot_wider(names_from=side, values_from=c(n,T)) %>%
  filter(!is.na(T_A), !is.na(T_B), n_A>=15, n_B>=15) %>%
  mutate(pair=paste(sub("Drosera_","",a), sub("Drosera_","",b), sep="/"),
         diff=round(T_A-T_B,3)) %>%
  select(pair, n_A, n_B, T_A, T_B, diff) %>% arrange(T_A)
print(as.data.frame(cmp), row.names=FALSE)
if (nrow(cmp) >= 4) {
  cat(sprintf("\n  Spearman T_A vs T_B: rho = %.3f\n",
              cor(cmp$T_A, cmp$T_B, method="spearman")))
  cat(sprintf("  median |T_A - T_B| = %.3f   vs median T_A = %.3f\n",
              median(abs(cmp$diff)), median(cmp$T_A)))
  cat(sprintf("  Wilcoxon paired: p = %.3g\n",
              suppressWarnings(wilcox.test(cmp$T_A, cmp$T_B, paired=TRUE)$p.value)))
}
cat("\n  T_A ~ T_B across pairs = the A and B donors split at the same times,\n")
cat("  i.e. the two lineages travelled together. That resolves the B-side\n")
cat("  polytomy without needing a separate B analysis.\n")
cat("  Systematic offset = the donors have different divergence histories.\n")

hr("4. implied ordering, both lineages")
cat("  All depths rate-corrected and on one scale.\n\n")
print(as.data.frame(bind_rows(
  cmp %>% transmute(event=paste0("speciation ", pair), A=T_A, B=T_B),
  res %>% filter(cls=="cross-AB", band==0.15) %>%
    transmute(event=paste0("A/B split (", sub("Drosera_","",sp), ")"), A=T, B=NA_real_)
)) %>% arrange(desc(coalesce(A,B))), row.names=FALSE)
write_csv(cmp, "trees/qc/A_vs_B_times.csv")
write_csv(res, "trees/qc/band_sweep.csv")
