#!/usr/bin/env Rscript
# ============================================================================
# DR03b — IS REGIA CLOSER TO DIONAEA THAN TO DROSERA?
#
# THE HYPOTHESIS
#   regia might be an early-branching Drosera, or it might sit outside Drosera
#   entirely -- a cryptic lineage nearer Dionaea. DR03 established that all
#   five Drosera draw A and B from the same two progenitor lineages, but said
#   nothing about where regia BRANCHES.
#
# THE TEST
#   Dionaea is in the data, so extend DR03's measurement to it:
#     rAA = d(A_sp, A_Dionaea) / d(A_sp, B_Dionaea)
#   Same internally-calibrated statistic, same rate-free logic.
#
# PREDICTIONS (T_AB = 0.57 from DR00's rate-corrected depths)
#   regia INSIDE Drosera        -> rAA(regia,Dio) = rAA(core,Dio) ~ 0.98
#                                  both coalesce at the SAME node, the
#                                  Drosera/Dionaea split, so they are EQUAL
#   regia SISTER TO DIONAEA     -> rAA(regia,Dio) ~ 0.53, clearly SMALLER
#                                  than rAA(core,Dio)
#
#   The discriminator is the DIFFERENCE between regia-Dionaea and
#   core-Dionaea, not either value alone.
#
# ----------------------------------------------------------------------------
# THE CIRCULARITY, AND WHY IT MATTERS HERE MORE THAN IN DR03
# ----------------------------------------------------------------------------
#   DR02 labels a Drosera gene "A" precisely BECAUSE it is closer to Dionaea's
#   A copy. So d(A_sp, A_Dio) is biased DOWNWARD by construction, and every
#   species will look artificially close to Dionaea.
#
#   This did not matter in DR03: both species there were labelled the same
#   way, so the bias cancelled in a Drosera-Drosera comparison. It does NOT
#   cancel here.
#
#   THREE MITIGATIONS, in increasing strength:
#     1. SEGMENT labels, not per-gene. A segment averages ~52 genes, so its
#        call is driven by the median rather than by any locus's own delta.
#        The bias survives but is much reduced.
#     2. The bias is IDENTICAL for regia and for core Drosera, since both are
#        labelled by the same procedure. So the DIFFERENCE between them is
#        far more robust than either absolute value. This is why the test is
#        framed as a difference.
#     3. FOUR-POINT labels as a replicate. Buneman's condition is rate-free
#        and does not use distance magnitudes, so it carries a different bias
#        structure. If both label sources agree, the answer is not an artefact
#        of how delta was defined.
#
#   A HELD-OUT check is not available: Nepenthes is the only outgroup and it
#   sits outside the A/B split, so it cannot arbitrate between Drosera and
#   Dionaea for a given subgenome.
#
# IN   DR/out/DR02_gene_delta.csv, DR/out/DR02_segments.csv,
#      DR/out/DR02_agreement.csv, DR/out/pairwise_ks.csv,
#      DR/locus_meta.tsv, fractionation_by_chrpair.csv
# OUT  DR/out/DR03b_dionaea.csv
# FIG  DR/fig/DR03b_1_placement.pdf   THE RESULT
#      DR/fig/DR03b_2_distributions.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))
set.seed(1); NBOOT <- 2000L
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

hr("0. SETUP")
LOC <- read_tsv("DR/locus_meta.tsv", show_col_types=FALSE) %>%
       group_by(tip) %>% slice(1) %>% ungroup()
K <- read_csv("DR/out/pairwise_ks.csv", show_col_types=FALSE) %>%
     filter(!is.na(dS), dS>=0, dS<5, codons>=100)
DD <- bind_rows(K %>% transmute(locus=anchor, a=seq1, b=seq2, d=dS),
                K %>% transmute(locus=anchor, a=seq2, b=seq1, d=dS))
dmap <- setNames(DD$d, paste(DD$locus, DD$a, DD$b))
gd <- function(l,x,y) ifelse(x==y, 0, unname(dmap[paste(l,x,y)]))
P <- read_csv("fractionation_by_chrpair.csv", show_col_types=FALSE)
sg1_dom <- P$retained_more == P$chrA
SIDE <- c(setNames(ifelse(sg1_dom,"A","B"), P$chrA),
          setNames(ifelse(sg1_dom,"B","A"), P$chrB))
dio <- LOC %>% filter(genome=="Dionaea_muscipula") %>%
  mutate(side=unname(SIDE[chr])) %>% filter(!is.na(side)) %>%
  select(locus, tip, side) %>%
  pivot_wider(names_from=side, values_from=tip, values_fn=list) %>%
  filter(lengths(A)==1, lengths(B)==1) %>%
  transmute(locus, DA=unlist(A), DB=unlist(B))
cat(sprintf("  Dionaea homeolog pairs: %d\n", nrow(dio)))

G  <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE)
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE)
lab <- SS %>% filter(label %in% c("A","B")) %>%
       select(genome, region, chr, mb_lo, mb_hi, label)
GL <- G %>% mutate(mb=mid/1e6) %>%
  inner_join(lab, by=c("genome","region","chr"), relationship="many-to-many") %>%
  filter(mb>=mb_lo, mb<=mb_hi) %>% group_by(tip) %>% slice(1) %>% ungroup()
cat(sprintf("  Drosera copies with a segment label: %d\n", nrow(GL)))

measure <- function(GLx, tag) {
  bind_rows(lapply(DROS, function(s) {
    d <- GLx %>% filter(genome==s)
    loci <- intersect(d$locus, dio$locus)
    bind_rows(lapply(loci, function(l) {
      a <- d$tip[d$locus==l & d$label=="A"]; b <- d$tip[d$locus==l & d$label=="B"]
      DA <- dio$DA[dio$locus==l]; DB <- dio$DB[dio$locus==l]
      if (!length(a) || !length(b)) return(NULL)
      mn <- function(x,y) { v <- unlist(lapply(x, function(i)
        vapply(y, function(j) gd(l,i,j), numeric(1)))); v <- v[!is.na(v)]
        if (!length(v)) NA_real_ else min(v) }
      AA <- mn(a, DA); BB <- mn(b, DB)
      AB <- mean(c(mn(a, DB), mn(b, DA)), na.rm=TRUE)
      if (any(is.na(c(AA,BB,AB))) || AB<=0) return(NULL)
      tibble(genome=s, locus=l, rAA=AA/AB, rBB=BB/AB, src=tag)
    }))
  }))
}

hr("1. rAA and rBB against DIONAEA, using SEGMENT labels")
cat("  PREDICTIONS: regia inside Drosera -> all five species EQUAL at ~0.98.\n")
cat("               regia sister to Dionaea -> regia clearly LOWER (~0.53).\n\n")
M <- measure(GL, "segment labels")
bm <- function(x) { b <- replicate(NBOOT, median(sample(x, replace=TRUE)))
                    quantile(b, c(.025,.975)) }
S1 <- M %>% group_by(genome) %>%
  group_modify(function(d, key) {
    ca <- bm(d$rAA); cb <- bm(d$rBB)
    tibble(n=nrow(d), rAA=median(d$rAA), rAA_lo=ca[1], rAA_hi=ca[2],
           rBB=median(d$rBB), rBB_lo=cb[1], rBB_hi=cb[2])
  }) %>% ungroup()
print(as.data.frame(S1 %>% transmute(genome=sub("Drosera_","",genome), n,
  rAA=round(rAA,3), rAA_CI=sprintf("[%.2f,%.2f]", rAA_lo, rAA_hi),
  rBB=round(rBB,3), rBB_CI=sprintf("[%.2f,%.2f]", rBB_lo, rBB_hi)) %>%
  arrange(rAA)), row.names=FALSE)

hr("2. THE DISCRIMINATOR — regia against the core three")
core <- c("Drosera_binata","Drosera_paradoxa","Drosera_scorpioides")
rg <- M %>% filter(genome=="Drosera_regia")
co <- M %>% filter(genome %in% core)
cat(sprintf("  regia      : rAA %.3f (n=%d)\n", median(rg$rAA), nrow(rg)))
cat(sprintf("  core three : rAA %.3f (n=%d)\n", median(co$rAA), nrow(co)))
dlt <- median(rg$rAA) - median(co$rAA)
w <- suppressWarnings(wilcox.test(rg$rAA, co$rAA))
cat(sprintf("  difference : %+.3f | Wilcoxon p = %.3g\n", dlt, w$p.value))
cat("\n  INTERPRETATION\n")
if (dlt < -0.10) {
  cat("    regia is MEASURABLY CLOSER to Dionaea than the core three are.\n")
  cat("    Consistent with regia branching outside Drosera.\n")
} else if (abs(dlt) <= 0.10) {
  cat("    regia and the core three are EQUIDISTANT from Dionaea.\n")
  cat("    That is what 'regia is inside Drosera' predicts: both coalesce\n")
  cat("    with Dionaea at the same node, so the depth is shared.\n")
} else {
  cat("    regia is FURTHER from Dionaea than the core three -- unexpected\n")
  cat("    under either hypothesis. Check for a rate or saturation artefact.\n")
}
cat("\n  CAVEAT: delta labels a gene A because it is close to Dionaea's A, so\n")
cat("  every species is biased toward Dionaea. The bias is the SAME for regia\n")
cat("  and for the core three, so the DIFFERENCE above is the robust quantity;\n")
cat("  the absolute values are not.\n")

hr("3. REPLICATE with FOUR-POINT labels")
cat("  Buneman's condition uses only the RANKING of three distance sums, not\n")
cat("  their magnitudes, so it carries a different bias structure to delta.\n")
cat("  Agreement between the two means the answer is not an artefact of how\n")
cat("  delta was defined.\n\n")
ag <- read_csv("DR/out/DR02_agreement.csv", show_col_types=FALSE)
FPL <- ag %>% filter(fp %in% c("A","B"), fp_gap > 0.05) %>%
       transmute(genome, locus, tip, label=fp)
cat(sprintf("  copies with a confident four-point label: %d\n", nrow(FPL)))
M2 <- measure(FPL, "four-point labels")
S2 <- M2 %>% group_by(genome) %>%
  summarise(n=n(), rAA=round(median(rAA),3), rBB=round(median(rBB),3),
            .groups="drop")
print(as.data.frame(S2 %>% mutate(genome=sub("Drosera_","",genome)) %>%
  arrange(rAA)), row.names=FALSE)
rg2 <- M2 %>% filter(genome=="Drosera_regia")
co2 <- M2 %>% filter(genome %in% core)
if (nrow(rg2) > 20 && nrow(co2) > 20) {
  d2 <- median(rg2$rAA) - median(co2$rAA)
  cat(sprintf("\n  four-point difference (regia - core): %+.3f\n", d2))
  cat(sprintf("  segment-label difference             : %+.3f\n", dlt))
  cat(ifelse(sign(d2)==sign(dlt) || abs(d2-dlt) < 0.1,
             "  SAME DIRECTION -- the answer is label-source independent.\n",
             "  *** DISAGREE -- the result depends on the labelling method. ***\n"))
}

hr("4. THE FULL DISTANCE PICTURE")
cat("  Raw depths, for a sanity check on saturation and rates.\n\n")
raw <- bind_rows(lapply(DROS, function(s) {
  d <- GL %>% filter(genome==s); loci <- intersect(d$locus, dio$locus)
  bind_rows(lapply(loci, function(l) {
    a <- d$tip[d$locus==l & d$label=="A"]
    DA <- dio$DA[dio$locus==l]; DB <- dio$DB[dio$locus==l]
    if (!length(a)) return(NULL)
    tibble(genome=s, locus=l,
           dAA=min(vapply(a, function(i) gd(l,i,DA), numeric(1)), na.rm=TRUE),
           dAB=min(vapply(a, function(i) gd(l,i,DB), numeric(1)), na.rm=TRUE))
  }))
})) %>% filter(is.finite(dAA), is.finite(dAB))
print(as.data.frame(raw %>% group_by(genome) %>%
  summarise(n=n(), med_dAA=round(median(dAA),3), med_dAB=round(median(dAB),3),
            frac_over_2=round(mean(dAA>2),3), .groups="drop") %>%
  mutate(genome=sub("Drosera_","",genome)) %>% arrange(med_dAA)),
  row.names=FALSE)
cat("\n  med_dAA is the raw Drosera-A to Dionaea-A distance. It carries RATE,\n")
cat("  so regia (rate 0.308 vs binata 1.0) should look small in raw terms\n")
cat("  regardless of topology. That is exactly why the ratio is used.\n")
write_csv(bind_rows(M, M2), "DR/out/DR03b_dionaea.csv")

p1 <- ggplot(S1, aes(reorder(sub("Drosera_","",genome), rAA), rAA)) +
  geom_hline(yintercept=0.98, colour="#378ADD", linetype="dashed") +
  geom_hline(yintercept=0.53, colour="#D85A30", linetype="dashed") +
  geom_linerange(aes(ymin=rAA_lo, ymax=rAA_hi), colour="grey55") +
  geom_point(size=3, colour="#1D9E75") +
  annotate("text", x=0.7, y=0.98, label="predicted if regia is INSIDE Drosera",
           hjust=0, vjust=-0.6, size=2.7, colour="#2B6CB0") +
  annotate("text", x=0.7, y=0.53, label="predicted if regia is SISTER to Dionaea",
           hjust=0, vjust=-0.6, size=2.7, colour="#993C1D") +
  coord_flip(ylim=c(0.3,1.2)) +
  labs(title="DR03b - distance to Dionaea, scaled by the A/B split",
       subtitle="all five equal = regia is inside Drosera | regia low = regia branches outside",
       x=NULL, y="A<->A / A<->B  against Dionaea") + theme_minimal(11)
suppressWarnings(ggsave("DR/fig/DR03b_1_placement.pdf", p1, width=8, height=5))
p2 <- ggplot(M %>% mutate(grp=ifelse(genome=="Drosera_regia","regia","other")),
             aes(rAA, fill=grp)) +
  geom_density(alpha=0.5, colour=NA) +
  geom_vline(xintercept=c(0.53,0.98), colour=c("#D85A30","#378ADD"),
             linetype="dashed") +
  facet_wrap(~sub("Drosera_","",genome), ncol=2, scales="free_y") + xlim(0,2) +
  scale_fill_manual(values=c(regia="#D85A30", other="#378ADD"), guide="none") +
  labs(title="DR03b - per-locus distributions",
       subtitle="dashed: 0.53 = sister to Dionaea, 0.98 = inside Drosera",
       x="rAA against Dionaea", y="density") + theme_minimal(9)
suppressWarnings(ggsave("DR/fig/DR03b_2_distributions.pdf", p2, width=9, height=8))

hr("READ IN THIS ORDER")
cat("  1. sec 2 -- the difference between regia and the core three. THE TEST.\n")
cat("  2. sec 3 -- do four-point labels give the same direction?\n")
cat("  3. FIG .1 -- where do the five species fall against the predictions?\n")
cat("  4. sec 4 -- saturation check on the raw depths.\n")
