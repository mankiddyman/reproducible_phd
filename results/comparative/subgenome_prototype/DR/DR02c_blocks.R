#!/usr/bin/env Rscript
# ============================================================================
# DR02c — PROPAGATION BY BLOCK MEMBERSHIP
#
# WHY BOUNDING BOXES FAILED
#   mb_lo/mb_hi is min/max of a segment's voting genes. But an ancestral region
#   can appear in TWO separated blocks on one chromosome:
#     binata chr4_hap1:  chr7(64) chr10(16) chr4(78) ... chr7(19) chr3(68)
#   chr7's box therefore runs from gene 1 to gene 192 and ENCLOSES the chr10
#   and chr4 blocks between its two appearances. The boxes overlap; the actual
#   gene blocks do not.
#
# WHAT THE DATA SUPPORTS
#   Run-length encoding of region along each chromosome: median run 45 genes,
#   only 4% single-gene runs, 45 chromosomes carrying more than one region.
#   Regions form long coherent blocks. So block membership is well defined and
#   is the right unit -- not a coordinate range.
#
# THE METHOD, with no invented parameters
#   1. Genes WITH a region (from gene_delta) are ordered along the chromosome.
#   2. rle() gives the blocks. Boundaries sit midway between the last gene of
#      one run and the first of the next. Nothing is chosen; the data draws it.
#   3. Every gene in combBed falls in exactly one block -> gets that region.
#   4. The block inherits the label of the segments ITS OWN voting genes
#      belong to. If those disagree, the block is dropped, not guessed.
#
# WHAT THIS DOES NOT FIX
#   The 19/19 result conditions on regions having exactly 3 segments, which is
#   close to conditioning on the answer. 19 of 40 (species x region) pairs
#   qualify; capensis contributes none. That is a separate problem and is NOT
#   addressed here.
#
# IN   DR/out/DR02_segments.csv, DR/out/DR02_gene_delta.csv,
#      ../genespace/results/combBed.txt
# OUT  DR/out/DR02_propagated_blocks.csv
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
GSD <- file.path(dirname(getwd()), "genespace", "results")
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

hr("1. INPUTS")
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE) %>%
      filter(label %in% c("A","B"))
GD <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE) %>%
      mutate(mb = mid/1e6)
OLD <- read_csv("DR/out/DR02_propagated.csv", show_col_types=FALSE)
allg <- read.table(file.path(GSD,"combBed.txt"), header=TRUE, sep="\t",
                   quote="", comment.char="", stringsAsFactors=FALSE) %>%
  mutate(isRep=as.logical(isArrayRep)) %>% filter(is.na(isRep)|isRep) %>%
  filter(genome %in% DROS) %>%
  transmute(genome, gene=id, chr, mb=(start+end)/2e6)
cat(sprintf("  segments: %d | voting genes: %d | genes to label: %d\n",
            nrow(SS), nrow(GD), nrow(allg)))

hr("2. ARE SEGMENTS WITHIN ONE (genome, region, chr) NON-OVERLAPPING?")
cat("  They should be -- they come from changepoint splitting of one track.\n")
cat("  If they overlap, assigning a voting gene to a segment is ambiguous.\n\n")
bad <- 0
for (kk in unique(paste(SS$genome, SS$region, SS$chr))) {
  z <- SS[paste(SS$genome,SS$region,SS$chr)==kk, ]
  if (nrow(z) < 2) next
  z <- z[order(z$mb_lo), ]
  if (any(z$mb_lo[-1] < z$mb_hi[-nrow(z)])) bad <- bad + 1
}
cat(sprintf("  tracks with overlapping segments: %d\n", bad))

hr("3. ASSIGN EACH VOTING GENE TO ITS SEGMENT")
GD$seg <- NA_character_; GD$seglab <- NA_character_
for (i in seq_len(nrow(SS))) {
  ix <- which(GD$genome==SS$genome[i] & GD$chr==SS$chr[i] &
              GD$region==SS$region[i] &
              GD$mb >= SS$mb_lo[i] & GD$mb <= SS$mb_hi[i] & is.na(GD$seg))
  GD$seg[ix] <- SS$segment[i]; GD$seglab[ix] <- SS$label[i]
}
V <- GD %>% filter(!is.na(seg))
cat(sprintf("  voting genes on a called segment: %d of %d\n", nrow(V), nrow(GD)))

hr("4. BLOCKS FROM RUN-LENGTH ENCODING")
BLK <- do.call(rbind, lapply(split(V, paste(V$genome, V$chr)), function(z){
  z <- z[order(z$mb), ]
  r <- rle(paste(z$region, z$seglab))   # label change ends a block
  ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  do.call(rbind, lapply(seq_along(r$lengths), function(j){
    w <- z[starts[j]:ends[j], ]
    tl <- table(w$seglab)
    data.frame(genome=z$genome[1], chr=z$chr[1], region=w$region[1],
               blk=j, nvote=nrow(w),
               first_mb=min(w$mb), last_mb=max(w$mb),
               label=names(tl)[which.max(tl)],
               purity=max(tl)/sum(tl), stringsAsFactors=FALSE)
  }))
}))
cat(sprintf("  blocks: %d | median voting genes per block: %.0f\n",
            nrow(BLK), median(BLK$nvote)))
cat(sprintf("  blocks whose voting genes DISAGREE on label (<80%% pure): %d -- dropped\n",
            sum(BLK$purity < 0.8)))
BLK <- BLK[BLK$purity >= 0.8 & BLK$nvote >= 3, ]
cat(sprintf("  blocks kept: %d (%d A, %d B)\n", nrow(BLK),
            sum(BLK$label=="A"), sum(BLK$label=="B")))

hr("5. BLOCK BOUNDARIES — midway between neighbouring blocks")
BLK <- BLK[order(BLK$genome, BLK$chr, BLK$first_mb), ]
BLK$lo <- NA_real_; BLK$hi <- NA_real_
for (kk in unique(paste(BLK$genome, BLK$chr))) {
  ix <- which(paste(BLK$genome,BLK$chr)==kk)
  z <- BLK[ix, ]
  lo <- c(0, (z$last_mb[-nrow(z)] + z$first_mb[-1])/2)
  hi <- c((z$last_mb[-nrow(z)] + z$first_mb[-1])/2, Inf)
  BLK$lo[ix] <- lo; BLK$hi[ix] <- hi
}
cat("  boundaries drawn from the data; no threshold chosen.\n")
o <- 0
for (kk in unique(paste(BLK$genome,BLK$chr))) {
  z <- BLK[paste(BLK$genome,BLK$chr)==kk, ]
  if (nrow(z) > 1 && any(z$lo[-1] < z$hi[-nrow(z)] - 1e-9)) o <- o + 1
}
cat(sprintf("  chromosomes with overlapping BLOCKS: %d  (boxes gave 25)\n", o))

hr("6. PROPAGATE")
NEW <- allg %>%
  inner_join(BLK %>% select(genome, chr, region, label, lo, hi),
             by=c("genome","chr"), relationship="many-to-many") %>%
  filter(mb >= lo, mb < hi) %>%
  group_by(genome, gene) %>%
  summarise(chr=chr[1], mb=mb[1], region=region[1], label=label[1],
            nmatch=n(), .groups="drop")
cat(sprintf("  genes matching more than one block: %d (should be 0)\n",
            sum(NEW$nmatch > 1)))
NEW <- NEW %>% select(-nmatch) %>% mutate(mid = mb*1e6)
cat(sprintf("  genes assigned: %d\n", nrow(NEW)))

hr("7. OLD vs NEW")
cmp <- allg %>% count(genome, name="total") %>%
  left_join(OLD %>% count(genome, name="old"), by="genome") %>%
  left_join(NEW %>% count(genome, name="new"), by="genome") %>%
  mutate(old_pct=round(100*old/total,1), new_pct=round(100*new/total,1))
print(as.data.frame(cmp), row.names=FALSE)
j <- OLD %>% select(genome, gene, ol=label, or=region) %>%
     inner_join(NEW %>% select(genome, gene, nl=label, nr=region),
                by=c("genome","gene"))
cat(sprintf("\n  genes in both: %d | label flips: %d (%.2f%%) | region changes: %d (%.2f%%)\n",
            nrow(j), sum(j$ol!=j$nl), 100*mean(j$ol!=j$nl),
            sum(j$or!=j$nr), 100*mean(j$or!=j$nr)))

hr("8. THE TEST CASE")
for (nm in c("old","new")) {
  d <- if (nm=="old") OLD else NEW
  cat(sprintf("\n  %s -- binata chr5_hap1:\n", nm))
  print(as.data.frame(d %>% filter(genome=="Drosera_binata", chr=="chr5_hap1") %>%
        count(region, label)), row.names=FALSE)
}
write_csv(NEW, "DR/out/DR02_propagated_blocks.csv")
cat("\n  wrote DR/out/DR02_propagated_blocks.csv\n")
