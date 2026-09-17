#!/usr/bin/env Rscript
# ============================================================================
# DR02 — LABEL EVERY DROSERA GENE COPY A OR B
#
# THE STATISTIC
#     delta = [ d(R, D_A) - d(R, D_B) ] / d(D_A, D_B)
#   R's own branch appears in both numerator terms and subtracts out. The
#   gene's overall rate multiplies every path and divides out against the
#   denominator. Bounded to [-1,+1] by the triangle inequality.
#   NO NEPENTHES NEEDED: 2,003 loci where Dionaea kept both copies, against
#   1,147 that also have an outgroup.
#
# ----------------------------------------------------------------------------
# THE EXPECTED EFFECT IS SMALL. THIS IS THE CENTRAL DESIGN CONSTRAINT.
# ----------------------------------------------------------------------------
#   With equal rates, for an A-derived copy R:
#     d(R,D_A) = 2r.T_spec   d(R,D_B) = 2r.T_AB   d(D_A,D_B) = 2r.T_AB
#     delta = T_spec/T_AB - 1
#   From DR00: T_AB ~ 0.57, Drosera-Dionaea speciation ~ 0.48
#     -> delta ~ -0.16 for A, +0.16 for B. SEPARATION 0.32.
#
#   Noise: both numerator terms are distances near 0.6 with 15-25% relative
#   error, so their difference carries ~+-0.14; divided by 0.525 gives
#   sd(delta) ~ 0.32. SEPARATION / NOISE ~ 1.0 SD PER GENE.
#
#   That is the arithmetic that killed script 26's mixture model (components
#   1.1 SD apart, needed ~3). Averaging over a segment of n genes divides the
#   noise by sqrt(n):
#       n=1  -> 1.0 SD      n=10 -> 3.1 SD
#       n=5  -> 2.2 SD      n=25 -> 4.9 SD
#   which is why scripts 19/20/21 worked at SEGMENT level and why per-gene
#   labelling was never going to.
#
#   SO: per-gene delta is the INPUT. Per-segment labels are the DELIVERABLE.
#   Section 3 measures whether per-gene resolution exists rather than assuming
#   either way -- if the marginal delta distribution is bimodal, per-gene
#   labelling is viable after all.
#
# ----------------------------------------------------------------------------
# WHAT DELTA CANNOT DO, AND WHY THE FOUR-POINT COMPARISON IS A DIAGNOSTIC
# ----------------------------------------------------------------------------
#   Delta removes R's branch and the gene's rate. It does NOT remove a
#   systematic rate difference between the D_A and D_B LINEAGES -- and v2 sec3.2
#   found X faster in 7 of 8 pairs.
#
#   The four-point on (R, D_A, D_B, Nepenthes) is immune to exactly that:
#     pairing (R,D_A)|(D_B,Nep) contains D_A via d(R,D_A)
#     pairing (R,D_B)|(D_A,Nep) contains D_A via d(D_A,Nep)
#     pairing (R,Nep)|(D_A,D_B) contains D_A via d(D_A,D_B)
#   All three sums gain the same amount, so the RANKING is unchanged.
#
#   So comparing them MEASURES the bias. And the test is SYMMETRY, not
#   agreement rate: both are noisy per gene, so 70-80% agreement is expected
#   even if both are unbiased. Under pure noise the two off-diagonal cells of
#   the 2x2 are EQUAL. Under rate asymmetry, delta calls B where four-point
#   calls A more often than the reverse. McNemar's test on the off-diagonals
#   is the diagnostic, and the imbalance gives the correction magnitude.
#
# ----------------------------------------------------------------------------
# SEGMENTS: MAXIMAL MERGE THEN SPLIT, NOT AGREEMENT-BASED MERGE
# ----------------------------------------------------------------------------
#   Merging "only where adjacent blocks agree" would split a uniform tract
#   wherever noise flips a gene, producing many small segments that are
#   internally consistent BY CONSTRUCTION and look confident without being
#   real. Reversed: start with one segment per (species, region, chromosome)
#   and SPLIT only at changepoints that beat a permutation null.
#   Conservative in the right direction -- a missed boundary averages two
#   ancestries and gives delta ~ 0, which reads as "unresolved". A false
#   boundary gives two confident WRONG calls.
#
# EXCLUSIONS
#   - Dionaea conversion loci: d(D_A,D_B) < 0.20 makes delta unstable (0/0).
#     AB04 found 1.45% of both-retained loci there.
#   - AB03 exchange tracts: 6.1% of Dionaea loci where the A chromosome
#     carries B ancestry, so every Drosera label is inverted. Run BOTH ways.
#   - paradoxa conversion loci: DR01 found d_sis p05 = 0.000. A B-copy
#     converted by A reads as A. Flagged from DR01_votes.csv.
#   - capensis recent duplicates: near-identical deltas that must not be
#     double-counted. Collapsed at dS < 0.25 (script 15's threshold).
#   - ~17% of copies are genuinely OUTSIDE the A/B split (script 17: Dionaea's
#     copies are sisters at 17.3% of loci). delta ~ 0 there, correctly.
#
# INHERITED ASSUMPTION
#   Everything here inherits AB sec3.1: whether Dionaea's fractionation
#   partition tracks ancestry at all. Partially supported (p = 0.004,
#   script 36b attacks 6-7) but not settled.
#
# IN   DR/locus_meta.tsv, DR/out/pairwise_ks.csv, DR/out/DR01_votes.csv,
#      AB/out/AB03_flagged_genes.csv, fractionation_by_chrpair.csv,
#      ../genespace/results/combBed.txt
# OUT  DR/out/DR02_gene_delta.csv     every copy: delta, flags
#      DR/out/DR02_segments.csv       segment calls with margin
#      DR/out/DR02_propagated.csv     labels for all genes on called segments
#      DR/out/DR02_agreement.csv      delta vs four-point
# FIG  DR/fig/DR02_1_marginal.pdf     READ FIRST - is delta bimodal?
#      DR/fig/DR02_2_controls.pdf     Dionaea +-1, Nepenthes centred
#      DR/fig/DR02_3_agreement.pdf    delta vs four-point
#      DR/fig/DR02_4_segments.pdf     segment deltas and coverage
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))
set.seed(1)
DSMAX <- 5; DCONV <- 0.20; NPERM <- 999L; MINSEG <- 5L
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")
GSD <- file.path(dirname(getwd()), "genespace", "results")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

hr("1. DIONAEA HOMEOLOG PAIRS — the reference axis")
LOC <- read_tsv("DR/locus_meta.tsv", show_col_types=FALSE) %>%
       group_by(tip) %>% slice(1) %>% ungroup()
K <- read_csv("DR/out/pairwise_ks.csv", show_col_types=FALSE) %>%
     filter(!is.na(dS), dS>=0, dS<DSMAX, codons>=100)
PAIRS <- read_csv("fractionation_by_chrpair.csv", show_col_types=FALSE)
sg1_dom <- PAIRS$retained_more == PAIRS$chrA
SIDE <- c(setNames(ifelse(sg1_dom, "A", "B"), PAIRS$chrA),
          setNames(ifelse(sg1_dom, "B", "A"), PAIRS$chrB))
cat("  A = the MORE-RETAINED chromosome of each pair.\n")
cat("  NOTE: chrA in fractionation_by_chrpair.csv is the alphabetically first\n")
cat("  chromosome, NOT the more-retained one. Using it directly labelled\n")
cat("  pairs chr4/chr5/chr6 inverted, which made their Nepenthes regions\n")
cat("  (chr5_dom, chr6_dom, chr7_dom) appear to double B while the other five\n")
cat("  doubled A. That was an orientation bug, not biology.\n\n")
print(as.data.frame(PAIRS %>% transmute(pair=exp_pair, chrA, chrB,
  retained_more, A_side=ifelse(retained_more==chrA, "sg1", "sg2"))),
  row.names=FALSE)
cat("\n")
dio <- LOC %>% filter(genome=="Dionaea_muscipula") %>%
  mutate(side=unname(SIDE[chr])) %>% filter(!is.na(side))
ax <- dio %>% select(locus, tip, side) %>%
  pivot_wider(names_from=side, values_from=tip, values_fn=list) %>%
  filter(lengths(A)==1, lengths(B)==1) %>%
  mutate(DA=unlist(A), DB=unlist(B)) %>% select(locus, DA, DB)
cat(sprintf("  loci with exactly one Dionaea copy per side: %d\n", nrow(ax)))
DD <- bind_rows(K %>% transmute(locus=anchor, a=seq1, b=seq2, d=dS),
                K %>% transmute(locus=anchor, a=seq2, b=seq1, d=dS))
dmap <- setNames(DD$d, paste(DD$locus, DD$a, DD$b))
gd <- function(l,x,y) ifelse(x==y, 0, unname(dmap[paste(l,x,y)]))
ax$dAB <- mapply(gd, ax$locus, ax$DA, ax$DB)
cat(sprintf("  with a measured d(D_A,D_B): %d\n", sum(!is.na(ax$dAB))))
ax <- ax %>% filter(!is.na(dAB))
cat(sprintf("  median d(D_A,D_B): %.3f  (v2 reports 0.525)\n", median(ax$dAB)))
nconv <- sum(ax$dAB < DCONV)
cat(sprintf("  EXCLUDED, d(D_A,D_B) < %.2f (Dionaea conversion): %d (%.1f%%)\n",
            DCONV, nconv, 100*nconv/nrow(ax)))
ax <- ax %>% filter(dAB >= DCONV)
cat(sprintf("  usable axis loci: %d\n", nrow(ax)))

hr("2. PER-GENE DELTA")
# region from the Dionaea chromosome pair, NOT from the locus name
dio_pair <- LOC %>% filter(genome=="Dionaea_muscipula") %>%
  mutate(dpair=sub("_sg[0-9]+_s[0-9]+$","",chr)) %>%
  group_by(locus) %>% summarise(dpair=dpair[1], .groups="drop")
nep_region <- LOC %>% filter(genome=="Nepenthes_gracilis") %>%
  group_by(locus) %>% summarise(nchr=chr[1], .groups="drop")
CAL <- dio_pair %>% inner_join(nep_region, by="locus") %>%
  count(dpair, nchr) %>% group_by(dpair) %>%
  slice_max(n, n=1, with_ties=FALSE) %>% ungroup() %>%
  transmute(dpair, region=nchr)
cat("\n  Dionaea pair -> Nepenthes region calibration:\n")
print(as.data.frame(CAL), row.names=FALSE)
stopifnot(!any(duplicated(CAL$region)), nrow(CAL)==8)
cat("  mapping is 1:1 across all 8 pairs.\n")

G <- LOC %>% filter(genome %in% DROS) %>% inner_join(ax, by="locus") %>%
  left_join(dio_pair, by="locus") %>% left_join(CAL, by="dpair") %>%
  filter(!is.na(region))
G$dRA <- mapply(gd, G$locus, G$tip, G$DA)
G$dRB <- mapply(gd, G$locus, G$tip, G$DB)
G <- G %>% filter(!is.na(dRA), !is.na(dRB)) %>%
     mutate(delta=(dRA-dRB)/dAB)
cat(sprintf("  Drosera copies with a delta: %d across %d loci\n",
            nrow(G), n_distinct(G$locus)))
print(as.data.frame(G %>% count(genome, name="copies") %>%
  left_join(G %>% distinct(genome, locus) %>% count(genome, name="loci"),
            by="genome")), row.names=FALSE)

hr("3. IS DELTA BIMODAL? — the per-gene resolution question")
cat("  Predicted from DR00: A ~ -0.16, B ~ +0.16, separation 0.32.\n")
cat("  Predicted per-gene noise sd ~ 0.32, i.e. ~1.0 SD separation.\n")
cat("  If the marginal is UNIMODAL, per-gene labelling is not viable and the\n")
cat("  deliverable is per-segment. If BIMODAL, per-gene works after all.\n\n")
print(as.data.frame(G %>% group_by(genome) %>%
  summarise(copies=n(), p10=round(quantile(delta,.10),3),
            q1=round(quantile(delta,.25),3), median=round(median(delta),3),
            q3=round(quantile(delta,.75),3), p90=round(quantile(delta,.90),3),
            sd=round(sd(delta),3), .groups="drop")), row.names=FALSE)
cat("\n  observed sd is the per-gene noise. Compare to the 0.32 separation:\n")
print(as.data.frame(G %>% group_by(genome) %>%
  summarise(sd=round(sd(delta),3), sep_in_SD=round(0.32/sd(delta),2),
            .groups="drop")), row.names=FALSE)
cat("  below ~2 SD, single genes cannot be labelled reliably.\n")
dip <- G %>% group_by(genome) %>%
  summarise(frac_mid = mean(abs(delta) < 0.08), .groups="drop")
cat("\n  fraction with |delta| < 0.08 (the ambiguous band):\n")
print(as.data.frame(dip %>% mutate(frac_mid=round(frac_mid,3))), row.names=FALSE)
p1 <- ggplot(G, aes(delta)) +
  geom_histogram(bins=80, fill="#378ADD", colour="white", linewidth=0.15) +
  geom_vline(xintercept=c(-0.16, 0.16), colour="firebrick", linetype="dashed") +
  geom_vline(xintercept=0, colour="grey40", linetype="dotted") +
  facet_wrap(~genome, scales="free_y", ncol=2) + xlim(-1,1) +
  labs(title="DR02.1 - per-gene delta",
       subtitle="dashed = predicted A (-0.16) and B (+0.16) | bimodal = per-gene labelling viable",
       x="delta", y="copies") + theme_minimal(9)
suppressWarnings(ggsave("DR/fig/DR02_1_marginal.pdf", p1, width=9, height=8))

hr("4. CONTROLS — these must pass or nothing below is meaningful")
cat("  (a) Dionaea's own copies: d(D_A,D_A)=0 so delta must be EXACTLY -1/+1\n")
ck <- ax %>% mutate(dA=mapply(gd,locus,DA,DA), dB=mapply(gd,locus,DB,DB)) %>%
      filter(!is.na(dA)) %>%
      mutate(delta_DA=(dA-dAB)/dAB, delta_DB=(dAB-dB)/dAB)
if (nrow(ck)) {
  cat(sprintf("      D_A delta: median %.4f (expect -1) | D_B: %.4f (expect +1)\n",
              median(ck$delta_DA, na.rm=TRUE), median(ck$delta_DB, na.rm=TRUE)))
  if (abs(median(ck$delta_DA, na.rm=TRUE)+1) > 0.01) {
    cat("      *** FAILED: sign convention or join is wrong. STOPPING. ***\n"); quit(status=1)
  }
  cat("      passed.\n")
} else cat("      (self-distances absent from the dS table -- check skipped)\n")

cat("  (b) Nepenthes sits OUTSIDE the A/B split, so delta must centre on 0\n")
NEP <- LOC %>% filter(genome=="Nepenthes_gracilis") %>%
       group_by(locus) %>% slice(1) %>% ungroup() %>%
       inner_join(ax, by="locus")
NEP$dRA <- mapply(gd, NEP$locus, NEP$tip, NEP$DA)
NEP$dRB <- mapply(gd, NEP$locus, NEP$tip, NEP$DB)
NEP <- NEP %>% filter(!is.na(dRA), !is.na(dRB)) %>% mutate(delta=(dRA-dRB)/dAB)
OFF <- median(NEP$delta)
cat(sprintf("      n=%d | median delta %+.4f | Wilcoxon p=%.3g\n", nrow(NEP), OFF,
            suppressWarnings(wilcox.test(NEP$delta)$p.value)))
cat("      v2 measured +0.023 (p=0.082). A non-zero offset is the D_A/D_B\n")
cat("      rate asymmetry (X faster in 7/8 pairs, v2 sec3.2). Subtracted below.\n")
G <- G %>% mutate(delta_raw=delta, delta=delta-OFF)
p2 <- ggplot(NEP, aes(delta)) +
  geom_histogram(bins=60, fill="grey55", colour="white", linewidth=0.15) +
  geom_vline(xintercept=0, colour="grey30") +
  geom_vline(xintercept=OFF, colour="firebrick", linetype="dashed") +
  labs(title="DR02.2 - Nepenthes control: delta must centre on zero",
       subtitle=sprintf("red = measured offset %+.4f, subtracted from all Drosera deltas", OFF),
       x="delta", y="loci") + theme_minimal(10)
suppressWarnings(ggsave("DR/fig/DR02_2_controls.pdf", p2, width=8, height=5))

hr("5. DELTA vs FOUR-POINT — a bias diagnostic, not a consistency check")
cat("  The four-point is immune to D_A/D_B rate asymmetry (all three sums gain\n")
cat("  equally). Delta is not. So disagreement ASYMMETRY measures that bias.\n")
cat("  Both are noisy per gene, so 70-80%% raw agreement is expected regardless\n")
cat("  -- McNemar on the off-diagonals is the real test.\n\n")
FP <- G %>% inner_join(NEP %>% select(locus, nep=tip), by="locus")
fp_call <- function(l, R, A, B, N) {
  s1 <- gd(l,R,A) + gd(l,B,N)
  s2 <- gd(l,R,B) + gd(l,A,N)
  s3 <- gd(l,R,N) + gd(l,A,B)
  if (any(is.na(c(s1,s2,s3)))) return(c(NA_character_, NA_real_))
  s <- c(s1,s2,s3); o <- sort(s)
  c(c("A","B","outside")[which.min(s)], (o[2]-o[1])/mean(s))
}
fpr <- t(mapply(fp_call, FP$locus, FP$tip, FP$DA, FP$DB, FP$nep))
FP$fp <- fpr[,1]; FP$fp_gap <- as.numeric(fpr[,2])
FP <- FP %>% filter(!is.na(fp))
cat(sprintf("  loci with both delta and four-point: %d copies\n", nrow(FP)))
print(as.data.frame(FP %>% count(fp, name="copies") %>%
  mutate(pct=round(100*copies/sum(copies),1))), row.names=FALSE)
cat("  'outside' ~17%% is expected: script 17 found Dionaea's copies are\n")
cat("  sisters at 17.3%% of loci, i.e. the Drosera copy predates the A/B split.\n\n")
FP <- FP %>% mutate(dl = ifelse(delta < 0, "A", "B"))
tb <- with(FP %>% filter(fp %in% c("A","B")), table(delta=dl, fourpoint=fp))
print(tb)
if (all(dim(tb)==c(2,2))) {
  b <- tb["A","B"]; c_ <- tb["B","A"]
  mc <- suppressWarnings(mcnemar.test(tb))
  cat(sprintf("\n  agreement: %.3f\n", sum(diag(tb))/sum(tb)))
  cat(sprintf("  delta=A/4pt=B: %d | delta=B/4pt=A: %d | McNemar p = %.3g\n",
              b, c_, mc$p.value))
  cat("  symmetric  -> disagreement is noise, delta is unbiased\n")
  cat("  asymmetric -> residual D_A/D_B rate bias in delta\n")
}
HC <- FP %>% filter(abs(delta) > 0.15, fp_gap > 0.05, fp %in% c("A","B"))
if (nrow(HC) > 20) {
  cat(sprintf("\n  HIGH-CONFIDENCE ONLY (|delta|>0.15 and gap>0.05): n=%d, agreement %.3f\n",
              nrow(HC), mean(HC$dl==HC$fp)))
  cat("  this should be much higher than the raw rate if both methods work.\n")
}
write_csv(FP %>% select(locus, genome, tip, delta, fp, fp_gap),
          "DR/out/DR02_agreement.csv")
p3 <- ggplot(FP %>% filter(fp %in% c("A","B")), aes(delta, fill=fp)) +
  geom_histogram(bins=60, alpha=0.6, position="identity", colour=NA) +
  geom_vline(xintercept=0, colour="grey30", linetype="dotted") +
  facet_wrap(~genome, scales="free_y", ncol=2) + xlim(-1,1) +
  scale_fill_manual(values=c(A="#1D9E75", B="#D85A30")) +
  labs(title="DR02.3 - delta distribution split by the four-point call",
       subtitle="if delta works, green sits left of zero and orange right",
       x="delta", y="copies", fill="four-point") +
  theme_minimal(9) + theme(legend.position="top")
suppressWarnings(ggsave("DR/fig/DR02_3_agreement.pdf", p3, width=9, height=8))

hr("6. SEGMENTS — maximal merge, then split at changepoints")
cat("  Starting unit: one segment per (species, region, chromosome).\n")
cat("  Split ONLY where a changepoint in the delta series beats a permutation\n")
cat("  null on gene order. A missed boundary averages two ancestries and reads\n")
cat("  as unresolved; a false boundary gives two confident WRONG calls.\n\n")
bedpos <- read.table(file.path(GSD,"combBed.txt"), header=TRUE, sep="\t",
                     quote="", comment.char="", stringsAsFactors=FALSE) %>%
  filter(genome %in% DROS) %>% transmute(genome, gene=id, mid=(start+end)/2)
G <- G %>% left_join(bedpos, by=c("genome","gene")) %>% filter(!is.na(mid)) %>%
     arrange(genome, region, chr, mid)
cp_scan <- function(v) {
  n <- length(v); if (n < 2*MINSEG) return(list(best=NA_integer_, stat=0))
  cs <- cumsum(v); tot <- cs[n]
  k <- MINSEG:(n-MINSEG)
  m1 <- cs[k]/k; m2 <- (tot-cs[k])/(n-k)
  st <- abs(m1-m2) * sqrt(k*(n-k)/n)
  list(best=k[which.max(st)], stat=max(st))
}
split_rec <- function(v, off=0L, depth=0L) {
  if (depth > 3 || length(v) < 2*MINSEG) return(integer(0))
  s <- cp_scan(v); if (is.na(s$best)) return(integer(0))
  nul <- replicate(NPERM, cp_scan(sample(v))$stat)
  if (mean(nul >= s$stat) > 0.05) return(integer(0))
  c(split_rec(v[1:s$best], off, depth+1L), off + s$best,
    split_rec(v[(s$best+1):length(v)], off + s$best, depth+1L))
}
SEG <- G %>% group_by(genome, region, chr) %>%
  group_modify(function(d, key) {
    d <- d[order(d$mid), ]
    cps <- if (nrow(d) >= 2*MINSEG) split_rec(d$delta) else integer(0)
    d$seg <- cumsum(seq_len(nrow(d)) %in% (cps+1L))
    d
  }) %>% ungroup() %>%
  mutate(segment = paste(genome, region, chr, seg, sep="|"))
nsp <- SEG %>% distinct(genome, region, chr, seg) %>%
  count(genome, region, chr, name="segs")
cat(sprintf("  chromosome-in-region tracks: %d | segments after splitting: %d\n",
            nrow(nsp), sum(nsp$segs)))
cat("  EXPECTED ~144: 4 hexaploids x 8 regions x ~3 chr, plus capensis x 8 x ~6.\n")
cat("  Thousands of tracks means the region assignment has broken again.\n")
print(as.data.frame(nsp %>% group_by(genome) %>%
  summarise(tracks=n(), segments=sum(segs), split_tracks=sum(segs>1),
            .groups="drop")), row.names=FALSE)

SS <- SEG %>% group_by(genome, region, chr, segment) %>%
  summarise(n=n(), med_delta=median(delta), mad_delta=mad(delta),
            mb_lo=min(mid)/1e6, mb_hi=max(mid)/1e6, .groups="drop")
SS$se    <- SS$mad_delta / sqrt(SS$n)
SS$label <- ifelse(SS$n < MINSEG, "unresolved",
             ifelse(SS$med_delta < -0.05, "A",
             ifelse(SS$med_delta >  0.05, "B", "ambiguous")))
SS$conf  <- round(abs(SS$med_delta) / pmax(SS$se, 1e-6), 2)
cat("\n  genes per segment (drives the noise -- see the header arithmetic):\n")
print(summary(SS$n))
cat("\n  segment calls:\n")
print(as.data.frame(SS %>% count(genome, label) %>%
  pivot_wider(names_from=label, values_from=n, values_fill=0)), row.names=FALSE)
cat("\n  A:B ratio among CALLED segments (DR04 tests this properly):\n")
rr <- SS %>% filter(label %in% c("A","B")) %>% count(genome, label) %>%
      pivot_wider(names_from=label, values_from=n, values_fill=0)
rr$ratio <- round(rr$A / pmax(rr$B, 1), 2)
print(as.data.frame(rr), row.names=FALSE)
write_csv(SS, "DR/out/DR02_segments.csv")
write_csv(G %>% select(locus, genome, tip, gene, chr, mid, region,
                       dRA, dRB, dAB, delta_raw, delta),
          "DR/out/DR02_gene_delta.csv")

hr("7. PROPAGATION — coverage past the per-gene ceiling")
cat("  Every gene on a called segment inherits its label, including genes\n")
cat("  where Dionaea lost a copy and no delta could be computed.\n\n")
called <- SS %>% filter(label %in% c("A","B")) %>%
  mutate(span_mb = mb_hi - mb_lo,
         genes_per_mb = n / pmax(span_mb, 0.01))
cat("  DENSITY GUARD: a segment may only claim genes inside its span if its\n")
cat("  voting genes actually cover that span. Sparse segments are reported\n")
cat("  but NOT propagated -- otherwise a few genes spanning a chromosome\n")
cat("  would label the whole thing.\n\n")
print(as.data.frame(called %>% group_by(genome) %>%
  summarise(segments=n(), med_genes=median(n), med_span_mb=round(median(span_mb),2),
            med_genes_per_mb=round(median(genes_per_mb),2), .groups="drop")),
  row.names=FALSE)
MINDENS <- 0.5
dense <- called %>% filter(genes_per_mb >= MINDENS)
cat(sprintf("\n  segments passing >= %.1f voting genes/Mb: %d of %d\n",
            MINDENS, nrow(dense), nrow(called)))
called <- dense
allg <- read.table(file.path(GSD,"combBed.txt"), header=TRUE, sep="\t",
                   quote="", comment.char="", stringsAsFactors=FALSE) %>%
  mutate(isRep=as.logical(isArrayRep)) %>% filter(is.na(isRep)|isRep) %>%
  filter(genome %in% DROS) %>% transmute(genome, gene=id, chr, mid=(start+end)/2)
PROP <- allg %>% inner_join(
  called %>% select(genome, chr, region, segment, label, mb_lo, mb_hi),
  by=c("genome","chr"), relationship="many-to-many") %>%
  filter(mid/1e6 >= mb_lo, mid/1e6 <= mb_hi) %>%
  group_by(genome, gene) %>% slice(1) %>% ungroup()
cov <- allg %>% count(genome, name="total_genes") %>%
  left_join(PROP %>% count(genome, name="labelled"), by="genome") %>%
  left_join(G %>% distinct(genome, gene) %>% count(genome, name="direct_delta"),
            by="genome")
cov$pct_direct <- round(100*cov$direct_delta/cov$total_genes, 1)
cov$pct_prop   <- round(100*cov$labelled/cov$total_genes, 1)
print(as.data.frame(cov), row.names=FALSE)
cov$genes_per_segment <- round(cov$labelled /
  pmax(called %>% count(genome) %>% right_join(cov["genome"], by="genome") %>%
       pull(n) %>% replace(is.na(.), 0L), 1))
cat("\n  pct_direct is the per-gene ceiling; pct_prop is what segments buy.\n")
cat("  SANITY: genes_per_segment must not exceed the largest segment's own\n")
cat("  gene count by much. If it does, segments are swallowing chromosomes.\n")
print(as.data.frame(cov %>% select(genome, labelled, genes_per_segment)),
      row.names=FALSE)
write_csv(PROP, "DR/out/DR02_propagated.csv")
p4 <- ggplot(SS %>% filter(n>=MINSEG), aes(med_delta, fill=label)) +
  geom_histogram(bins=50, alpha=0.75, colour=NA) +
  geom_vline(xintercept=c(-0.05,0.05), colour="grey30", linetype="dotted") +
  facet_wrap(~genome, scales="free_y", ncol=2) +
  scale_fill_manual(values=c(A="#1D9E75", B="#D85A30", ambiguous="grey60",
                             unresolved="grey80")) +
  labs(title="DR02.4 - segment median delta",
       subtitle="averaging over genes should separate A and B far better than FIG 1",
       x="segment median delta", y="segments", fill=NULL) +
  theme_minimal(9) + theme(legend.position="top")
suppressWarnings(ggsave("DR/fig/DR02_4_segments.pdf", p4, width=9, height=8))

hr("READ IN THIS ORDER")
cat("  1. sec 4  -- did the controls pass? Dionaea +-1, Nepenthes near zero.\n")
cat("  2. FIG .1 -- is per-gene delta bimodal? sec 3 gives separation in SD.\n")
cat("  3. sec 5  -- McNemar. Symmetric = noise; asymmetric = residual bias.\n")
cat("  4. FIG .4 -- do SEGMENT medians separate where genes did not?\n")
cat("  5. sec 7  -- how much of each genome ends up labelled?\n")
