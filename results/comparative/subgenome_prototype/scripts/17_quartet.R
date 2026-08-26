#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"; DIO <- "Dionaea_muscipula"

k0 <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS > 0, codons >= 100)

build <- function(CEN) {
  k <- k0 %>% filter(dS < CEN)
  DD <- bind_rows(k %>% transmute(anchor, ga=seq1, gb=seq2, d=dS),
                  k %>% transmute(anchor, ga=seq2, gb=seq1, d=dS))
  fp <- k %>% filter(sp1==sp2, sp1!=NEP) %>%
        transmute(anchor, sp=sp1, S1=seq1, S2=seq2, dSS=dS) %>%
        group_by(anchor, sp) %>% filter(n()==1) %>% ungroup()
  nep <- k %>% filter(xor(sp1==NEP, sp2==NEP)) %>%
         transmute(anchor, N=ifelse(sp1==NEP, seq1, seq2)) %>% distinct()
  tg <- bind_rows(k %>% transmute(anchor, g=seq1, tsp=sp1),
                  k %>% transmute(anchor, g=seq2, tsp=sp2)) %>%
        distinct() %>% filter(tsp != NEP)
  q <- fp %>% inner_join(nep, by="anchor") %>%
       inner_join(tg %>% rename(Tg=g), by="anchor") %>% filter(tsp != sp) %>%
       left_join(DD %>% rename(S1=ga, Tg=gb, d_S1T=d), by=c("anchor","S1","Tg")) %>%
       left_join(DD %>% rename(S2=ga, Tg=gb, d_S2T=d), by=c("anchor","S2","Tg")) %>%
       left_join(DD %>% rename(S1=ga, N=gb,  d_S1N=d), by=c("anchor","S1","N")) %>%
       left_join(DD %>% rename(S2=ga, N=gb,  d_S2N=d), by=c("anchor","S2","N")) %>%
       left_join(DD %>% rename(Tg=ga, N=gb,  d_TN=d),  by=c("anchor","Tg","N")) %>%
       filter(!is.na(d_S1T), !is.na(d_S2T), !is.na(d_S1N), !is.na(d_S2N), !is.na(d_TN))
  q %>% mutate(
    sA = dSS   + d_TN,      # (S1,S2 | T,N)  copies are SISTERS
    sB = d_S1T + d_S2N,     # (S1,T  | S2,N) T allies with S1
    sC = d_S2T + d_S1N,     # (S2,T  | S1,N) T allies with S2
    m1 = pmin(sA,sB,sC), m3 = pmax(sA,sB,sC), m2 = sA+sB+sC-m1-m3,
    support = (m2-m1)/((sA+sB+sC)/3),
    call = case_when(sA<=sB & sA<=sC ~ "sisters",
                     sB<=sA & sB<=sC ~ "with_S1", TRUE ~ "with_S2"),
    interleaved = call != "sisters", CEN = CEN)
}

hr("1. quartet resolution by focal species  (null = 1/3 sisters)")
cat("  Four-point condition on (S1, S2, T, Nepenthes). Rate-free for additive dS.\n")
cat("  sisters >> 1/3 : lineage-specific WGD.  << 1/3 : copies interleave with T.\n\n")
q <- build(2.0)
cat(sprintf("  quartets: %d   loci: %d\n\n", nrow(q), n_distinct(q$anchor)))
tab <- q %>% group_by(sp) %>%
  summarise(quartets=n(), loci=n_distinct(anchor),
            sisters=round(mean(call=="sisters"),3),
            interleaved=round(mean(interleaved),3),
            med_support=round(median(support),4),
            unresolved=round(mean(support < 0.05),3), .groups="drop") %>%
  arrange(sisters)
tab$p_vs_third <- vapply(seq_len(nrow(tab)), function(i) {
  d <- q %>% filter(sp==tab$sp[i])
  signif(binom.test(sum(d$call=="sisters"), nrow(d), 1/3)$p.value, 3)
}, numeric(1))
print(as.data.frame(tab), row.names=FALSE)

hr("2. is this a rapid radiation? support = gap between best and second-best")
cat("  Near-zero support across the board = the internal branch is ~0 and the\n")
cat("  quartet genuinely cannot be resolved. That is the right-panel scenario.\n\n")
print(as.data.frame(q %>% group_by(sp) %>%
  summarise(p25=round(quantile(support,.25),4), med=round(median(support),4),
            p75=round(quantile(support,.75),4), p95=round(quantile(support,.95),4),
            .groups="drop")), row.names=FALSE)

hr("3. Dionaea focal, split by partner species")
print(as.data.frame(q %>% filter(sp==DIO) %>% group_by(tsp) %>%
  summarise(quartets=n(), sisters=round(mean(call=="sisters"),3),
            interleaved=round(mean(interleaved),3),
            med_support=round(median(support),4), .groups="drop") %>%
  arrange(sisters)), row.names=FALSE)
cat("\n  If regia interleaves most, it is nested deepest inside Dionaea's split.\n")

hr("4. censoring sweep — does saturation drive the call?")
cat("  Sum A uses a within-species distance; B and C use cross-species ones.\n")
cat("  Saturation compresses the larger values, biasing toward 'interleaved'.\n\n")
sw <- bind_rows(lapply(c(1.2, 1.5, 2.0, 3.0), function(cc) {
  build(cc) %>% group_by(sp) %>%
    summarise(CEN=cc, quartets=n(), sisters=round(mean(call=="sisters"),3), .groups="drop")
}))
print(as.data.frame(sw %>% select(sp, CEN, quartets, sisters) %>%
  pivot_wider(names_from=CEN, values_from=c(quartets, sisters))), row.names=FALSE)
cat("\n  Stable across thresholds = real. Drifting = saturation.\n")

hr("5. conservative: one call per locus (quartets share S1,S2,N)")
pl <- q %>% group_by(sp, anchor) %>%
  summarise(n_q=n(), f_sis=mean(call=="sisters"), .groups="drop") %>%
  mutate(loc_call = ifelse(f_sis > 0.5, "sisters", "interleaved"))
res <- pl %>% group_by(sp) %>%
  summarise(loci=n(), sisters=round(mean(loc_call=="sisters"),3), .groups="drop") %>%
  arrange(sisters)
res$p <- vapply(seq_len(nrow(res)), function(i) {
  d <- pl %>% filter(sp==res$sp[i])
  signif(binom.test(sum(d$loc_call=="sisters"), nrow(d), 1/3)$p.value, 3)
}, numeric(1))
print(as.data.frame(res), row.names=FALSE)

q %>% mutate(sp=sub("Drosera_","D_",sub("Dionaea_muscipula","DIONAEA",sp))) %>%
  ggplot(aes(support)) + geom_histogram(bins=60, fill="#378ADD", colour="white") +
  geom_vline(xintercept=0.05, linetype="dashed", colour="firebrick") +
  facet_wrap(~sp, scales="free_y") + xlim(0, 0.6) +
  labs(title="FIG50 - quartet resolution strength",
       subtitle="mass piled at zero = unresolvable internal branch",
       x="support (gap between best and second-best sum)", y="quartets") +
  theme_minimal(10) -> p
suppressWarnings(ggsave("trees/qc/FIG50_quartet_support.pdf", p, width=10, height=6))
write_csv(q %>% select(anchor, sp, tsp, dSS, sA, sB, sC, support, call),
          "trees/qc/quartet_calls.csv")
cat("\nWROTE: trees/qc/quartet_calls.csv  FIG50_quartet_support.pdf\n")
