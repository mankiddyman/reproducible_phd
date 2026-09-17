#!/usr/bin/env Rscript
# ============================================================================
# DRtest — ASSERTIONS THE PHASING PRODUCT MUST SATISFY
#
# Written after a day in which the same class of bug appeared four times:
#   1. slice(1) gave every gene in overlapping bounding boxes to one segment
#   2. rle(region) merged an A tract and a B tract that shared a region
#   3. the riparian re-derived tracts instead of using the ones already computed
#   4. a "sanity check" asserted something the data structure made impossible
# Each was found by looking at a figure, not by any test. That is backwards.
#
# Run this after ANY change to DR02, DR02c, or a figure that reads them.
# Exit status is non-zero if any assertion fails.
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
FAIL <- 0
ok <- function(cond, msg, detail="") {
  if (isTRUE(cond)) cat(sprintf("  PASS  %s\n", msg))
  else { cat(sprintf("  FAIL  %s\n", msg)); if (nzchar(detail)) cat("        ", detail, "\n")
         FAIL <<- FAIL + 1 }
}
hr <- function(t) cat(sprintf("\n=== %s ===\n", t))

S  <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE)
GD <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE)
PB <- read_csv("DR/out/DR02_propagated_blocks.csv", show_col_types=FALSE)
AB <- S %>% filter(label %in% c("A","B"))

hr("1. every called segment is represented")
miss <- AB %>%
  anti_join(PB %>% distinct(genome, chr, region, label),
            by=c("genome","chr","region","label"))
ok(nrow(miss)==0, "no called (genome,chr,region,label) is missing from the product",
   paste(head(paste(miss$genome, miss$chr, miss$region, miss$label), 5), collapse=" | "))

hr("2. a chromosome carrying BOTH labels shows both")
both <- AB %>% count(genome, chr, label) %>%
  count(genome, chr, name="nlab") %>% filter(nlab > 1)
bad <- both %>% left_join(
  PB %>% count(genome, chr, label) %>% count(genome, chr, name="plab"),
  by=c("genome","chr")) %>% filter(is.na(plab) | plab < 2)
ok(nrow(bad)==0,
   sprintf("all %d chromosomes with A and B segments show both", nrow(both)),
   paste(head(paste(bad$genome, bad$chr), 5), collapse=" | "))

hr("3. the specific case that broke twice")
r5 <- PB %>% filter(genome=="Drosera_regia", chr=="chr5_collapsed")
ok(nrow(r5) > 0, "regia chr5_collapsed has propagated genes")
ok(n_distinct(r5$label) == 2,
   "regia chr5_collapsed carries BOTH A and B (segments: A 0-10.2Mb, B 10.3-19.0Mb)",
   paste(names(table(r5$label)), table(r5$label), collapse=" "))
b5 <- PB %>% filter(genome=="Drosera_binata", chr=="chr5_hap1")
ok(n_distinct(b5$region) >= 2,
   "binata chr5_hap1 carries both chr3_dom and chr6_dom",
   paste(unique(b5$region), collapse=" "))

hr("4. no gene assigned twice")
d <- PB %>% count(genome, gene) %>% filter(n > 1)
ok(nrow(d)==0, "each gene appears exactly once", sprintf("%d duplicates", nrow(d)))

hr("5. labels agree with the segment calls")
chk <- PB %>% inner_join(AB %>% select(genome, chr, region, seg_label=label),
                         by=c("genome","chr","region"),
                         relationship="many-to-many") %>%
  group_by(genome, chr, region) %>%
  summarise(prop_labs=paste(sort(unique(label)), collapse=""),
            seg_labs=paste(sort(unique(seg_label)), collapse=""), .groups="drop") %>%
  filter(prop_labs != seg_labs)
ok(nrow(chk)==0, "propagated labels match the segment calls per (genome,chr,region)",
   paste(head(sprintf("%s %s %s: prop=%s seg=%s", chk$genome, chk$chr,
                      chk$region, chk$prop_labs, chk$seg_labs), 5), collapse=" | "))

hr("6. coverage is plausible")
cov <- PB %>% count(genome, name="n")
ok(all(cov$n > 5000), "every species has >5000 labelled genes",
   paste(cov$genome, cov$n, collapse=" "))

hr("7. THE INVARIANT THAT WOULD HAVE CAUGHT THE ORIGINAL BUG")
cat("  Every analysis script joins genes to segments on (genome, region, chr).\n")
cat("  If a gene can match more than one segment under that key, slice(1)\n")
cat("  starts choosing arbitrarily and every downstream result is at risk.\n")
cat("  Dropping region from the key gives 24.9%% multi-matches -- that was\n")
cat("  the bug in DR02_label.R line 391.\n\n")
.lab <- AB %>% select(genome, region, chr, mb_lo, mb_hi, segment, label)
.mm <- GD %>% mutate(mb=mid/1e6) %>%
  inner_join(.lab, by=c("genome","region","chr"), relationship="many-to-many") %>%
  filter(mb>=mb_lo, mb<=mb_hi) %>% count(tip, name="nmatch")
ok(sum(.mm$nmatch > 1) == 0,
   "no tip matches >1 segment under (genome, region, chr)",
   sprintf("%d tips multi-match -- slice(1) is choosing arbitrarily",
           sum(.mm$nmatch > 1)))

hr("8. the constitution is untouched")
k <- table(paste(AB$genome, AB$region))
three <- names(k)[k==3]
res <- t(sapply(three, function(x){
  z <- AB[paste(AB$genome,AB$region)==x, ]
  c(A=sum(z$label=="A"), B=sum(z$label=="B")) }))
n21 <- sum(res[,"A"]==2 & res[,"B"]==1)
ok(length(three)==19 && n21==18,
   "19 three-segment regions, 18 at 2A:1B",
   sprintf("found %d regions, %d at 2A:1B", length(three), n21))

cat(sprintf("\n%s\n%d FAILURE%s\n%s\n", strrep("=",60), FAIL,
            ifelse(FAIL==1,"","S"), strrep("=",60)))
quit(status = ifelse(FAIL > 0, 1, 0))
