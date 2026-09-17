#!/usr/bin/env Rscript
# ============================================================================
# DR22 — IS THE COPY-NUMBER CEILING REAL, OR TRUNCATED BY OUR LABELLING?
#
# THE PROBLEM WITH DR21
#   A_p95 = 2 says A copies never exceed 2 at 95% of loci. But only 4.8-30% of
#   loci actually CARRY 2 A copies -- binata 422 of 4530, capensis 217 of 4498.
#   So the ceiling is inferred from a thin upper tail.
#
#   binata has 2950 loci with exactly ONE A copy. Either
#     (a) the second A copy was lost to fractionation   -> ceiling stands
#     (b) the second A copy EXISTS but was not labelled -> counts truncated,
#         and the whole distribution is an artefact of coverage
#
# THE TEST
#   locus_meta.tsv lists every tip, labelled or not. For a locus short of the
#   predicted copy number, count TIPS against LABELS:
#     tips == labels  -> the copy is genuinely absent from the annotation (a)
#     tips >  labels  -> the copy is present but we failed to label it (b)
#
#   If (b) is common, DR21's ceiling is an artefact and the constitution claim
#   cannot rest on it.
#
# IN   DR/locus_meta.tsv, DR/out/DR02_propagated_blocks.csv
# OUT  DR/out/DR22_truncation.csv
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")
PRED <- c(Drosera_binata=2, Drosera_paradoxa=2, Drosera_regia=2,
          Drosera_scorpioides=2, Drosera_capensis=4)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

LOC <- read_tsv("DR/locus_meta.tsv", show_col_types=FALSE) %>% filter(genome %in% DROS)
PB  <- read_csv("DR/out/DR02_propagated_blocks.csv", show_col_types=FALSE)
key <- setNames(PB$label, paste(PB$genome, PB$gene))
LOC$lab <- unname(key[paste(LOC$genome, LOC$gene)])

hr("1. TIPS vs LABELS, per locus per species")
T <- LOC %>% group_by(genome, locus) %>%
  summarise(tips=n(), labelled=sum(!is.na(lab)),
            nA=sum(lab=="A", na.rm=TRUE), nB=sum(lab=="B", na.rm=TRUE),
            .groups="drop") %>%
  mutate(sp=sub("Drosera_","",genome), unlab=tips-labelled)
cat(sprintf("  locus x species rows: %d\n", nrow(T)))
print(as.data.frame(T %>% group_by(sp) %>%
  summarise(loci=n(), median_tips=median(tips), median_labelled=median(labelled),
            pct_fully_labelled=round(100*mean(unlab==0),1), .groups="drop")),
  row.names=FALSE)

hr("2. THE KEY QUESTION")
cat("  For loci BELOW the predicted A count, is the shortfall explained by\n")
cat("  missing tips, or by tips we failed to label?\n\n")
T$predA <- unname(PRED[T$genome])
short <- T %>% filter(nA < predA)
res <- short %>% group_by(sp) %>%
  summarise(n=n(),
            no_spare_tip = sum(unlab == 0),
            has_spare_tip = sum(unlab > 0),
            pct_explained_by_absence = round(100*mean(unlab==0),1),
            .groups="drop")
print(as.data.frame(res), row.names=FALSE)
cat("\n  no_spare_tip  = every tip at that locus IS labelled, so the missing\n")
cat("                  A copy is genuinely not in the annotation -> real loss\n")
cat("  has_spare_tip = an unlabelled tip exists that COULD be the missing copy\n")
cat("                  -> the count may be truncated by our labelling\n")

hr("3. WHERE DO THE UNLABELLED TIPS SIT?")
cat("  If unlabelled tips are on chromosomes with no called segment, that is\n")
cat("  a coverage gap. If they are on labelled chromosomes, something else.\n\n")
U <- LOC %>% filter(is.na(lab))
cat(sprintf("  unlabelled tips: %d of %d (%.1f%%)\n",
            nrow(U), nrow(LOC), 100*nrow(U)/nrow(LOC)))
if (nrow(U)) {
  labelled_chr <- PB %>% distinct(genome, chr) %>% mutate(has=TRUE)
  U <- U %>% left_join(labelled_chr, by=c("genome","chr"))
  print(as.data.frame(U %>% group_by(sp=sub("Drosera_","",genome)) %>%
    summarise(unlabelled=n(),
              on_labelled_chr=sum(!is.na(has)),
              on_unlabelled_chr=sum(is.na(has)), .groups="drop")),
    row.names=FALSE)
  cat("\n  top chromosomes carrying unlabelled tips:\n")
  print(as.data.frame(U %>% count(genome, chr, sort=TRUE) %>% head(12)),
        row.names=FALSE)
}

hr("4. CORRECTED CEILING — assume every unlabelled tip is a missing copy")
cat("  Worst case: add all unlabelled tips to the A count and see whether the\n")
cat("  ceiling moves. If it stays at 2, DR21 holds regardless.\n\n")
T$nA_max <- T$nA + T$unlab
print(as.data.frame(T %>% group_by(sp) %>%
  summarise(A_p95 = quantile(nA,.95), A_p95_worst = quantile(nA_max,.95),
            A_p99 = quantile(nA,.99), A_p99_worst = quantile(nA_max,.99),
            .groups="drop")), row.names=FALSE)
cat("\n  If A_p95_worst is still 2 for the core four and 4 for capensis, the\n")
cat("  constitution holds even under the most pessimistic assumption.\n")
write_csv(T, "DR/out/DR22_truncation.csv")
