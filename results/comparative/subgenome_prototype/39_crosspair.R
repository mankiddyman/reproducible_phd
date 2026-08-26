#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
PAIR <- c(chr1_sg1_s5="p1", chr1_sg2_s9="p1", chr2_sg1_s3="p2", chr2_sg2_s11="p2",
          chr3_sg1_s7="p3", chr3_sg2_s15="p3", chr4_sg1_s2="p4", chr4_sg2_s12="p4",
          chr5_sg1_s10="p5", chr5_sg2_s16="p5", chr6_sg1_s6="p6", chr6_sg2_s8="p6",
          chr7_sg1_s1="p7", chr7_sg2_s13="p7", chr8_sg1_s4="p8", chr8_sg2_s14="p8")
X <- c("chr1_sg1_s5","chr2_sg1_s3","chr3_sg1_s7","chr4_sg2_s12",
       "chr5_sg2_s16","chr6_sg2_s8","chr7_sg1_s1","chr8_sg1_s4")

d <- read_csv("distances.csv", show_col_types=FALSE) %>%
  filter(!is.na(dS_dio1_dio2), dS_dio1_dio2 > 0) %>%
  mutate(P1=unname(PAIR[chr1]), P2=unname(PAIR[chr2]),
         cls = case_when(chr1==chr2 ~ "same chromosome",
                         P1==P2     ~ "within pair",
                         TRUE       ~ "CROSS PAIR"),
         xy1 = ifelse(chr1 %in% X,"X","Y"), xy2 = ifelse(chr2 %in% X,"X","Y"))
cat("NOTE: distances.csv uses NG86 dS (distances.py), which runs ~6% BELOW YN00.\n")
cat("Dionaea homeologs are 0.493 here vs 0.525 in ks/pairwise_ks.csv, so the\n")
cat("homeolog window is taken as 0.28-0.80 rather than 0.30-0.85.\n")
LO <- 0.28; HI <- 0.80

hr("1. counts by arrangement")
print(as.data.frame(d %>% count(cls, name="n") %>% mutate(pct=round(100*n/sum(n),1))),
      row.names=FALSE)
cat("\n  v2 sec3.1 reported 699 within-pair, 191 cross-pair, 16 same-chromosome.\n")

hr("2. THE TEST — is the cross-pair class at homeolog depth?")
cat("  Translocated homeologs keep their original divergence, so they should sit\n")
cat("  in the same window as within-pair homeologs. Ancient dispersed paralogs\n")
cat("  or orthology errors should not.\n\n")
print(as.data.frame(d %>% group_by(cls) %>%
  summarise(n=n(), p05=round(quantile(dS_dio1_dio2,.05),3),
            p25=round(quantile(dS_dio1_dio2,.25),3),
            med=round(median(dS_dio1_dio2),3),
            p75=round(quantile(dS_dio1_dio2,.75),3),
            frac_in_window=round(mean(dS_dio1_dio2>LO & dS_dio1_dio2<HI),3),
            .groups="drop")), row.names=FALSE)
cat(sprintf("\n  within-pair frac_in_window is the benchmark: those ARE homeologs.\n"))

hr("3. Nepenthes check — are both copies equally diverged from the outgroup?")
cat("  A translocated homeolog keeps its distance to Nepenthes. A spurious\n")
cat("  orthology call or an ancient paralog need not.\n\n")
print(as.data.frame(d %>% filter(!is.na(dS_nep_dio1), !is.na(dS_nep_dio2)) %>%
  mutate(asym = abs(dS_nep_dio1 - dS_nep_dio2)/pmax(dS_nep_dio1, dS_nep_dio2)) %>%
  group_by(cls) %>%
  summarise(n=n(), med_nep1=round(median(dS_nep_dio1),3),
            med_nep2=round(median(dS_nep_dio2),3),
            med_asym=round(median(asym),3), .groups="drop")), row.names=FALSE)
cat("\n  Similar med_asym across classes = cross-pair copies are proper homeologs.\n")

hr("4. which chromosome pairs exchange, and is it concentrated?")
xp <- d %>% filter(cls=="CROSS PAIR") %>%
  mutate(a=pmin(P1,P2), b=pmax(P1,P2)) %>% count(a,b, name="n") %>% arrange(desc(n))
cat(sprintf("  distinct pair-combinations: %d of %d possible\n\n", nrow(xp), 28))
print(as.data.frame(head(xp, 15)), row.names=FALSE)
cat("\n  Concentrated in a few combinations = real translocations.\n")
cat("  Spread across all 28 = orthology noise.\n")

hr("5. do cross-pair copies respect the X/Y partition?")
cat("  If a homeolog moved to another chromosome pair but kept its ancestry,\n")
cat("  the two copies should still be one X and one Y.\n\n")
print(as.data.frame(d %>% filter(cls=="CROSS PAIR") %>%
  count(xy1, xy2, name="n") %>%
  mutate(type=ifelse(xy1==xy2, "same side", "OPPOSITE sides"))), row.names=FALSE)
cat("\n  Mostly opposite sides = ancestry preserved through the translocation.\n")
cat("  50/50 = the X/Y labels say nothing about these copies.\n")

hr("6. cross-pair cases at homeolog depth — the candidate list")
cand <- d %>% filter(cls=="CROSS PAIR", dS_dio1_dio2>LO, dS_dio1_dio2<HI) %>%
  transmute(nep, chr1, chr2, xy1, xy2, dS=round(dS_dio1_dio2,3)) %>% arrange(dS)
cat(sprintf("  candidates: %d\n\n", nrow(cand)))
if (nrow(cand)) {
  print(as.data.frame(head(cand, 25)), row.names=FALSE)
  cat(sprintf("\n  of these, opposite X/Y sides: %d (%.0f%%)\n",
              sum(cand$xy1!=cand$xy2), 100*mean(cand$xy1!=cand$xy2)))
}
write_csv(cand, "trees/qc/crosspair_candidates.csv")

p <- d %>% ggplot(aes(dS_dio1_dio2)) +
  geom_histogram(bins=50, fill="#378ADD", colour="white") +
  geom_vline(xintercept=c(LO,0.493,HI), linetype="dashed", colour="firebrick") +
  facet_wrap(~cls, scales="free_y", ncol=1) + xlim(0,3) +
  labs(title="FIG54 - Dionaea copy pairs by chromosomal arrangement (NG86 dS)",
       subtitle="dashed: 0.28 / 0.493 homeolog depth / 0.80", x="dS", y="pairs") +
  theme_minimal(10)
suppressWarnings(ggsave("trees/qc/FIG54_arrangement_ng86.pdf", p, width=8, height=8))
cat("\nWROTE: trees/qc/crosspair_candidates.csv  FIG54_arrangement_ng86.pdf\n")
