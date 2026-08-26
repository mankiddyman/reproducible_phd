#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); NPERM <- 999
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"; DIO <- "Dionaea_muscipula"
CORE <- c("Drosera_binata","Drosera_paradoxa","Drosera_scorpioides")
TH <- 0.3

k <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
     filter(!is.na(dS), dS>0, dS<5, codons>=100)

## Dionaea loci with exactly one homeolog pair
dio <- k %>% filter(sp1==DIO, sp2==DIO) %>% transmute(anchor, D1=seq1, D2=seq2, dDD=dS) %>%
       group_by(anchor) %>% filter(n()==1) %>% ungroup() %>% filter(dDD > 0.05)
cat(sprintf("Dionaea 2-copy loci usable: %d\n", nrow(dio)))

cx <- k %>% filter(xor(sp1==DIO, sp2==DIO), sp1!=NEP, sp2!=NEP) %>%
  transmute(anchor,
            dg = ifelse(sp1==DIO, seq1, seq2),
            og = ifelse(sp1==DIO, seq2, seq1),
            sp = ifelse(sp1==DIO, sp2, sp1), dS)

D <- dio %>%
  inner_join(cx, by=c("anchor","D1"="dg")) %>% rename(dA=dS) %>%
  inner_join(cx, by=c("anchor","D2"="dg","og","sp")) %>% rename(dB=dS) %>%
  mutate(delta = (dA - dB)/dDD) %>%
  filter(is.finite(delta), abs(delta) <= 1.001)
cat(sprintf("scored copies: %d across %d loci, %d species\n",
            nrow(D), n_distinct(D$anchor), n_distinct(D$sp)))

hr("1. Delta distribution per species")
cat("  Mass at both extremes = carries both Dionaea subgenome affinities.\n")
cat("  Mass at 0            = sits OUTSIDE the Dionaea A/B split.\n")
cat("  Mass at one extreme  = carries one lineage only.\n\n")
print(as.data.frame(D %>% group_by(sp) %>%
  summarise(copies=n(), loci=n_distinct(anchor),
            p10=round(quantile(delta,.10),2), med=round(median(delta),2),
            p90=round(quantile(delta,.90),2),
            neg=round(mean(delta < -TH),3), mid=round(mean(abs(delta) <= TH),3),
            pos=round(mean(delta > TH),3),
            med_abs=round(median(abs(delta)),3), .groups="drop") %>%
  arrange(desc(med_abs))), row.names=FALSE)
cat(sprintf("\n  med_abs near 1 = copies place cleanly on one side or the other.\n"))
cat("  med_abs near 0 = the species is unplaceable w.r.t. Dionaea's split.\n")

hr("2. do a species' OWN copies split across the two Dionaea subgenomes?")
cat("  Works at any copy number — this is what unblocks regia (hexaploid, 3 subgenomes).\n\n")
loc <- D %>% group_by(sp, anchor) %>%
  summarise(kc=n(), lo=min(delta), hi=max(delta),
            split = any(delta < -TH) & any(delta > TH),
            spread = max(delta)-min(delta), .groups="drop") %>% filter(kc >= 2)

set.seed(1)
obs <- loc %>% group_by(sp) %>%
  summarise(loci=n(), mean_k=round(mean(kc),2),
            f_split=round(mean(split),3), med_spread=round(median(spread),3), .groups="drop")
nul <- bind_rows(lapply(split(D, D$sp), function(d) {
  ka <- d %>% count(anchor) %>% filter(n>=2)
  if (nrow(ka) < 20) return(NULL)
  st <- vapply(seq_len(NPERM), function(i) {
    dd <- d; dd$delta <- sample(dd$delta)
    z <- dd %>% group_by(anchor) %>% filter(n()>=2) %>%
         summarise(s = any(delta < -TH) & any(delta > TH),
                   sp2 = max(delta)-min(delta), .groups="drop")
    c(mean(z$s), median(z$sp2))
  }, numeric(2))
  tibble(sp=d$sp[1], null_split=round(mean(st[1,]),3), null_spread=round(mean(st[2,]),3),
         hi_split=round(quantile(st[1,],.975),3))
}))
res <- obs %>% left_join(nul, by="sp") %>%
  mutate(verdict = case_when(is.na(null_split) ~ "n/a",
                             f_split > hi_split ~ "copies SPLIT",
                             TRUE ~ "not above null")) %>%
  arrange(desc(f_split))
print(as.data.frame(res), row.names=FALSE)
cat("\n  Permutation shuffles delta within species, keeping the marginal\n")
cat("  distribution and destroying only the within-locus pairing.\n")

hr("3. regia by copy number")
print(as.data.frame(loc %>% filter(sp=="Drosera_regia") %>%
  group_by(kc) %>% summarise(loci=n(), f_split=round(mean(split),3),
                             med_spread=round(median(spread),3), .groups="drop")),
  row.names=FALSE)
cat("\n  If 3-copy regia loci split at a high rate, the hexaploid carries\n")
cat("  more than one ancestral lineage and the 2x2 test was simply the wrong shape.\n")

hr("4. core-trio calibration")
print(as.data.frame(res %>% filter(sp %in% CORE) %>% select(sp, loci, f_split, null_split)),
      row.names=FALSE)

D %>% mutate(sp = sub("Drosera_","D_",sp)) %>%
  ggplot(aes(delta)) + geom_histogram(binwidth=0.08, fill="#7F77DD", colour="white") +
  geom_vline(xintercept=c(-TH,0,TH), linetype="dashed", colour="grey40") +
  facet_wrap(~sp, scales="free_y") +
  labs(title="FIG49 - placement of each Drosera copy on the Dionaea homeolog pair",
       subtitle="-1 allies with D1, 0 outside the split, +1 allies with D2",
       x="delta", y="copies") + theme_minimal(10) -> p
ggsave("trees/qc/FIG49_delta.pdf", p, width=10, height=6)
write_csv(D, "trees/qc/delta_placement.csv"); write_csv(loc, "trees/qc/delta_by_locus.csv")
cat("\nWROTE: trees/qc/{delta_placement,delta_by_locus}.csv  FIG49_delta.pdf\n")
