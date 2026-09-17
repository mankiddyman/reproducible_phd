#!/usr/bin/env Rscript
# ============================================================================
# DR08b — WHY DOES THE CONTROL GIVE 0.27 AND NOT 0.173?
#
# Script 17 asked whether Dionaea's A and B copies are sisters. Headline
# 17.3% against a 1/3 null, p = 4.1e-19 -- the evidence for allopolyploidy.
# DR08 re-ran that quartet as a control and got 0.27.
#
# BUT script 17 also swept dS censoring from 1.2 to 3.0 and reported the
# statistic STABLE at 0.228-0.254. DR08 used dS < 5, looser than any point in
# that sweep, and got 0.27 -- just above the top of the range. So censoring
# alone may explain it.
#
# Four candidate causes, tested one at a time:
#   1. dS CENSORING     DR08 used 5; script 17 swept 1.2-3.0
#   2. GENE SET         old wgd7 (861 loci) vs new DR set (990)
#   3. THIRD TIP        one Drosera copy vs every Drosera copy
#   4. GAP FILTER       DR08 drops near-ties at gap <= 0.02, which RAISES the
#                       sister rate if near-ties are disproportionately
#                       non-sister
#
# WHY THIS MATTERS: the same quartet engine produced DR05e's 0.523 (regia's
# placement) and DR07's 0.000-0.008. A systematic bias would propagate.
#
# OUT  DR/out/DR08b_reconcile.csv
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

LOC <- read_tsv("DR/locus_meta.tsv", show_col_types=FALSE) %>%
       group_by(tip) %>% slice(1) %>% ungroup()
P <- read_csv("fractionation_by_chrpair.csv", show_col_types=FALSE)
s1 <- P$retained_more == P$chrA
SIDE <- c(setNames(ifelse(s1,"A","B"), P$chrA), setNames(ifelse(s1,"B","A"), P$chrB))
dio <- LOC %>% filter(genome=="Dionaea_muscipula") %>%
  mutate(side=unname(SIDE[chr])) %>% filter(!is.na(side)) %>%
  select(locus, tip, side) %>%
  pivot_wider(names_from=side, values_from=tip, values_fn=list) %>%
  filter(lengths(A)==1, lengths(B)==1) %>%
  transmute(locus, DioA=unlist(A), DioB=unlist(B))
nep <- LOC %>% filter(genome=="Nepenthes_gracilis") %>%
       group_by(locus) %>% slice(1) %>% ungroup() %>% transmute(locus, Nep=tip)
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")
drt <- LOC %>% filter(genome %in% DROS) %>%
       transmute(locus, tip, sp=sub("Drosera_","",genome))
base <- dio %>% inner_join(nep, by="locus")

run_at <- function(dsmax, gapcut, third, ks_file) {
  K <- read_csv(ks_file, show_col_types=FALSE) %>%
       filter(!is.na(dS), dS>0, dS<dsmax, codons>=100)
  DD <- bind_rows(K %>% transmute(locus=anchor, a=seq1, b=seq2, d=dS),
                  K %>% transmute(locus=anchor, a=seq2, b=seq1, d=dS))
  dm <- setNames(DD$d, paste(DD$locus, DD$a, DD$b))
  g <- function(l,x,y) unname(dm[paste(l,x,y)])
  rows <- bind_rows(lapply(base$locus, function(l) {
    A <- base$DioA[base$locus==l]; B <- base$DioB[base$locus==l]
    N <- base$Nep[base$locus==l]
    T3 <- if (third=="binata") drt$tip[drt$locus==l & drt$sp=="binata"]
          else drt$tip[drt$locus==l]
    if (!length(T3)) return(NULL)
    bind_rows(lapply(T3, function(t3) {
      S <- c(g(l,A,B)+g(l,t3,N), g(l,A,t3)+g(l,B,N), g(l,A,N)+g(l,B,t3))
      if (any(is.na(S))) return(NULL)
      o <- sort(S)
      tibble(locus=l, sis=which.min(S)==1, gap=(o[2]-o[1])/mean(S))
    }))
  }))
  if (!nrow(rows)) return(NULL)
  r <- rows %>% filter(gap > gapcut)
  tibble(dS_max=dsmax, gap_cut=gapcut, third_tip=third,
         set=ifelse(grepl("DR/", ks_file), "new DR", "old wgd7"),
         n=nrow(r), sisters=round(mean(r$sis),3))
}

hr("1. dS CENSORING — script 17 swept 1.2 to 3.0 and got 0.228-0.254")
res <- bind_rows(lapply(c(1.2,1.5,2.0,2.5,3.0,5.0), function(d)
  run_at(d, 0.02, "binata", "DR/out/pairwise_ks.csv")))
print(as.data.frame(res), row.names=FALSE)
cat("\n  If the sister rate FALLS toward 0.23-0.25 as dS_max tightens, the\n")
cat("  discrepancy is censoring and the code is fine.\n")

hr("2. GAP FILTER — do near-ties bias the rate?")
res2 <- bind_rows(lapply(c(0, 0.01, 0.02, 0.05), function(gc)
  run_at(3.0, gc, "binata", "DR/out/pairwise_ks.csv")))
print(as.data.frame(res2), row.names=FALSE)
cat("\n  A rising sister rate with a stricter gap cut means near-ties are\n")
cat("  disproportionately NON-sister, so filtering them inflates the estimate.\n")

hr("3. THIRD TIP — one Drosera copy vs all of them")
res3 <- bind_rows(
  run_at(3.0, 0.02, "binata", "DR/out/pairwise_ks.csv"),
  run_at(3.0, 0.02, "all",    "DR/out/pairwise_ks.csv"))
print(as.data.frame(res3), row.names=FALSE)
cat("\n  'all' gives many quartets per locus, which are NOT independent --\n")
cat("  script 17 may have done this. It changes the estimate and inflates n.\n")

hr("4. GENE SET — old wgd7 vs the new DR build")
if (file.exists("ks/pairwise_ks.csv")) {
  res4 <- bind_rows(
    run_at(3.0, 0.02, "binata", "ks/pairwise_ks.csv"),
    run_at(3.0, 0.02, "binata", "DR/out/pairwise_ks.csv"))
  print(as.data.frame(res4), row.names=FALSE)
} else cat("  old ks/pairwise_ks.csv not found\n")

hr("VERDICT")
best <- res %>% filter(dS_max <= 3.0)
cat(sprintf("  sister rate at script 17's censoring range: %s\n",
            paste(sprintf("%.3f", best$sisters), collapse=", ")))
cat("  script 17 reported 0.228-0.254 across that same sweep.\n\n")
cat("  OVERLAPPING -> the code is sound. DR08's 0.27 came from the looser\n")
cat("    dS < 5, and its core-species numbers should be re-read at dS < 3.\n")
cat("  NOT OVERLAPPING -> a real implementation difference, and DR05e/DR07\n")
cat("    need re-checking too, since they share this engine.\n")
write_csv(bind_rows(res, res2, res3), "DR/out/DR08b_reconcile.csv")
