#!/usr/bin/env Rscript
# ============================================================================
# DR02b — PROPAGATION, REDONE ON TRACTS INSTEAD OF BOUNDING BOXES
#
# THE BUG IN DR02_label.R SECTION 7
#   PROP <- allg %>% inner_join(called, by=c("genome","chr")) %>%
#           filter(mid/1e6 >= mb_lo, mid/1e6 <= mb_hi) %>%
#           group_by(genome, gene) %>% slice(1)
#
#   mb_lo/mb_hi are min/max of a segment's VOTING genes -- a bounding box, not
#   a contiguous tract. A segment with a few scattered genes gets a box across
#   the whole chromosome. On binata chr5_hap1:
#       chr3_dom  0.07-47.74 Mb   (79 voting genes)
#       chr6_dom  5.83-44.97 Mb   (50 voting genes, entirely INSIDE the above)
#   Every gene in 5.8-45 Mb falls in both boxes. slice(1) gives them all to
#   whichever row sorted first: chr3_dom took 865 genes, chr6_dom got ZERO.
#   20 of 180 called segments end up with no genes at all.
#
# WHY IT DOES NOT INVALIDATE THE PHASING
#   Segment CALLS come from the voting genes, before this step. 19/19 regions
#   at 2:1, 18/19 doubling A, 0.97 cross-species agreement -- all count
#   SEGMENTS. Every downstream script (DR03, DR05, DR07, DR08, DR09, DR11)
#   reads DR02_segments.csv and does its own (genome, region, chr) lookup;
#   none reads this file. Verified by grep.
#
# THE FIX
#   Split each segment at gaps in its voting genes, so a tract is a real
#   contiguous stretch rather than a bounding box. Overlap between tracts then
#   becomes rare and the tie-break stops mattering.
#
#   THE GAP THRESHOLD IS THE ONE REMAINING JUDGEMENT CALL, so it is swept
#   rather than asserted, and the sensitivity is printed. Genes are assigned
#   only INSIDE a tract -- no nearest-neighbour extrapolation, no distance
#   cutoff, no rule invented to make a number look better.
#
# IN   DR/out/DR02_segments.csv, DR/out/DR02_gene_delta.csv,
#      ../genespace/results/combBed.txt
# OUT  DR/out/DR02_propagated_v2.csv   (original left untouched)
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
GSD <- file.path(dirname(getwd()), "genespace", "results")
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")
MINDENS <- 0.5
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

hr("1. INPUTS")
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE)
GD <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE)
OLD <- read_csv("DR/out/DR02_propagated.csv", show_col_types=FALSE)
called <- SS %>% filter(label %in% c("A","B")) %>%
  mutate(span_mb=mb_hi-mb_lo, genes_per_mb=n/pmax(span_mb,0.01)) %>%
  filter(genes_per_mb >= MINDENS)
cat(sprintf("  called segments passing the density guard: %d\n", nrow(called)))
allg <- read.table(file.path(GSD,"combBed.txt"), header=TRUE, sep="\t",
                   quote="", comment.char="", stringsAsFactors=FALSE) %>%
  mutate(isRep=as.logical(isArrayRep)) %>% filter(is.na(isRep)|isRep) %>%
  filter(genome %in% DROS) %>%
  transmute(genome, gene=id, chr, mb=(start+end)/2e6)
cat(sprintf("  genes to label: %d\n", nrow(allg)))

hr("2. HOW MUCH DO BOUNDING BOXES OVERLAP?")
ovl <- function(tab, lo, hi) {
  n <- 0
  for (g in unique(tab$genome)) for (cc in unique(tab$chr[tab$genome==g])) {
    z <- tab[tab$genome==g & tab$chr==cc, ]
    if (nrow(z) < 2) next
    z <- z[order(z[[lo]]), ]
    if (any(z[[lo]][-1] < z[[hi]][-nrow(z)])) n <- n + 1
  }
  n
}
cat(sprintf("  chromosomes where BOXES overlap: %d\n",
            ovl(as.data.frame(called), "mb_lo", "mb_hi")))

hr("3. GAP SWEEP — the one judgement call, made visible")
vote <- GD %>% inner_join(called %>% select(genome, chr, region, segment, label),
                          by=c("genome","chr","region"))
mk_tracts <- function(GAP) {
  vote %>% mutate(mb=mid/1e6) %>%
    group_by(segment) %>% arrange(mb, .by_group=TRUE) %>%
    mutate(piece = cumsum(c(0, diff(mb) > GAP))) %>%
    group_by(genome, chr, region, segment, label, piece) %>%
    summarise(t_lo=min(mb), t_hi=max(mb), nvote=n(), .groups="drop") %>%
    filter(nvote >= 3)
}
sweep <- bind_rows(lapply(c(0.25,0.5,1,2,5), function(GAP) {
  TR <- mk_tracts(GAP)
  P <- allg %>% inner_join(TR, by=c("genome","chr"),
                           relationship="many-to-many") %>%
       filter(mb >= t_lo, mb <= t_hi)
  amb <- P %>% group_by(genome, gene) %>%
         summarise(nlab=n_distinct(label), .groups="drop")
  tibble(gap=GAP, tracts=nrow(TR),
         chrom_overlap=ovl(as.data.frame(TR), "t_lo", "t_hi"),
         genes=n_distinct(paste(P$genome,P$gene)),
         pct=round(100*n_distinct(paste(P$genome,P$gene))/nrow(allg),1),
         ambiguous=sum(amb$nlab > 1),
         pct_ambig=round(100*mean(amb$nlab > 1),2))
}))
print(as.data.frame(sweep), row.names=FALSE)
cat("\n  ambiguous = genes falling in two tracts with DIFFERENT labels.\n")
cat("  Those are the only ones a tie-break would decide, so a low count\n")
cat("  means the choice of gap barely matters.\n")

hr("4. ASSIGN, AT gap = 1 Mb")
GAP <- 1
TR <- mk_tracts(GAP)
cat(sprintf("  tracts: %d (%d A, %d B)\n", nrow(TR),
            sum(TR$label=="A"), sum(TR$label=="B")))
P <- allg %>% inner_join(TR, by=c("genome","chr"),
                         relationship="many-to-many") %>%
     filter(mb >= t_lo, mb <= t_hi)
amb <- P %>% group_by(genome, gene) %>%
       summarise(nlab=n_distinct(label), .groups="drop") %>%
       filter(nlab > 1)
cat(sprintf("  genes in >1 tract with conflicting labels: %d -- DROPPED, not guessed\n",
            nrow(amb)))
NEW <- P %>% anti_join(amb, by=c("genome","gene")) %>%
  group_by(genome, gene) %>% slice_min(t_hi - t_lo, n=1, with_ties=FALSE) %>%
  ungroup() %>%
  transmute(genome, gene, chr, mid=mb*1e6, region, segment, label,
            mb_lo=t_lo, mb_hi=t_hi)
cat(sprintf("  genes assigned: %d\n", nrow(NEW)))

hr("5. OLD vs NEW")
cmp <- allg %>% count(genome, name="total") %>%
  left_join(OLD %>% count(genome, name="old"), by="genome") %>%
  left_join(NEW %>% count(genome, name="new"), by="genome") %>%
  mutate(old_pct=round(100*old/total,1), new_pct=round(100*new/total,1))
print(as.data.frame(cmp), row.names=FALSE)
cat(sprintf("\n  segments with ZERO genes -- old: %d | new: %d (of %d called)\n",
  nrow(anti_join(called, distinct(OLD, segment), by="segment")),
  nrow(anti_join(called, distinct(NEW, segment), by="segment")), nrow(called)))
j <- OLD %>% select(genome, gene, ol=label, or=region) %>%
     inner_join(NEW %>% select(genome, gene, nl=label, nr=region),
                by=c("genome","gene"))
cat(sprintf("  genes in both: %d | label flips: %d (%.2f%%) | region changes: %d (%.2f%%)\n",
            nrow(j), sum(j$ol!=j$nl), 100*mean(j$ol!=j$nl),
            sum(j$or!=j$nr), 100*mean(j$or!=j$nr)))

hr("6. THE TEST CASE — binata chr5_hap1")
for (nm in c("old","new")) {
  d <- if (nm=="old") OLD else NEW
  cat(sprintf("\n  %s:\n", nm))
  print(as.data.frame(d %>% filter(genome=="Drosera_binata", chr=="chr5_hap1") %>%
        count(region, label)), row.names=FALSE)
}

hr("7. DOES THE CONSTITUTION HOLD?")
cat("  Segment counts are untouched by this step -- shown for reassurance.\n\n")
k <- table(paste(SS$genome[SS$label %in% c("A","B")],
                 SS$region[SS$label %in% c("A","B")]))
three <- names(k)[k==3]
res <- t(sapply(three, function(x) {
  z <- SS[paste(SS$genome,SS$region)==x & SS$label %in% c("A","B"), ]
  c(A=sum(z$label=="A"), B=sum(z$label=="B")) }))
print(table(paste0(res[,"A"],"A:",res[,"B"],"B")))
write_csv(NEW, "DR/out/DR02_propagated_v2.csv")
cat("\n  wrote DR/out/DR02_propagated_v2.csv (original untouched)\n")
