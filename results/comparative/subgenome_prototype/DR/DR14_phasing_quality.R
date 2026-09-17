#!/usr/bin/env Rscript
# ============================================================================
# DR14 — HOW GOOD ARE THE PHASED SEGMENTS?
#
# The phasing reaches 82-94% coverage, but only 8-19% of genes have a directly
# measured delta. The rest inherit their segment's label. Three questions:
#
#   1. IS THE SIGNAL CLEAN, OR MIXED? A segment called A by 30 votes to 2 is
#      not the same as one called A by 18 to 14. The label is identical; the
#      confidence is not. This asks how many segments are near-unanimous and
#      how many are marginal.
#
#   2. HOW FAR IS A PROPAGATED GENE FROM THE NEAREST MEASURED ONE? A gene 50 kb
#      from a voting gene is well supported. One 20 Mb away is an extrapolation.
#      That distance is the honest per-gene confidence and belongs in the figure
#      as transparency rather than uniform colour.
#
#   3. DO THE SEGMENT BOUNDARIES FALL ON SYNTENY BOUNDARIES? The ancestral
#      REGION comes from GENESPACE (Dionaea chromosome pair -> Nepenthes _dom,
#      calibrated 1:1). But the splits WITHIN a region were drawn by a
#      changepoint scan on the delta series -- 179 tracks became 207 segments.
#      If those 28 splits land on GENESPACE block edges, that is independent
#      corroboration. If they fall mid-block, the segmentation is following
#      delta noise.
#
# IN   DR/out/DR02_gene_delta.csv, DR02_segments.csv, DR02_propagated.csv,
#      ../genespace/results/{combBed.txt,syntenicBlock_coordinates.csv}
# OUT  DR/out/DR14_segment_quality.csv, DR14_gene_confidence.csv
# FIG  DR/fig/DR14_1_vote_purity.pdf, DR14_2_distance.pdf, DR14_3_boundaries.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))
GSD <- file.path(dirname(getwd()), "genespace", "results")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

G  <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE)
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE)
cat(sprintf("  genes with a measured delta: %d | segments: %d\n", nrow(G), nrow(SS)))

hr("1. VOTE PURITY — is the signal clean or mixed?")
cat("  For each segment, how many of its genes agree with the call?\n")
cat("  Ambiguous genes (|delta| < 0.05) are counted separately -- they are\n")
cat("  neither agreement nor conflict, they are absence of signal.\n\n")
lab <- SS %>% select(genome, region, chr, mb_lo, mb_hi, segment, label, med_delta, n)
GA <- G %>% mutate(mb=mid/1e6) %>%
  inner_join(lab, by=c("genome","region","chr"), relationship="many-to-many") %>%
  filter(mb>=mb_lo, mb<=mb_hi) %>% group_by(tip) %>% slice(1) %>% ungroup() %>%
  mutate(gene_call = ifelse(abs(delta) < 0.05, "ambig",
                     ifelse(delta < 0, "A", "B")))
PUR <- GA %>% group_by(segment, genome, label) %>%
  summarise(votes=n(),
            agree = sum(gene_call == label),
            conflict = sum(gene_call %in% c("A","B") & gene_call != label),
            ambig = sum(gene_call == "ambig"), .groups="drop")
PUR$purity <- with(PUR, agree / pmax(agree + conflict, 1))
PUR$frac_ambig <- PUR$ambig / PUR$votes
PUR$p_binom <- mapply(function(a,c) if (a+c < 3) NA_real_ else
                        binom.test(a, a+c, 0.5)$p.value, PUR$agree, PUR$conflict)
cat("  purity = agreeing genes / (agreeing + conflicting), ambiguous excluded\n\n")
print(summary(PUR$purity))
cat("\n  by species:\n")
print(as.data.frame(PUR %>% group_by(genome) %>%
  summarise(segments=n(), median_votes=median(votes),
            median_purity=round(median(purity),3),
            median_frac_ambig=round(median(frac_ambig),3),
            near_unanimous=sum(purity > 0.8),
            marginal=sum(purity < 0.6),
            signif_p05=sum(p_binom < 0.05, na.rm=TRUE), .groups="drop")),
  row.names=FALSE)
cat("\n  MARGINAL segments (purity < 0.6) are the ones to worry about --\n")
cat("  their label is a coin flip dressed as a call.\n\n")
print(as.data.frame(PUR %>% filter(purity < 0.6) %>%
  transmute(segment, genome, label, votes, agree, conflict, ambig,
            purity=round(purity,3)) %>% arrange(purity)), row.names=FALSE)
write_csv(PUR, "DR/out/DR14_segment_quality.csv")

p1 <- ggplot(PUR, aes(purity, fill=genome)) +
  geom_histogram(bins=30, alpha=0.75, colour=NA) +
  geom_vline(xintercept=0.5, linetype="dashed", colour="firebrick") +
  facet_wrap(~genome, scales="free_y", ncol=2) +
  labs(title="DR14.1 - vote purity per segment",
       subtitle="0.5 = the genes disagree as often as they agree | ambiguous genes excluded",
       x="fraction of voting genes agreeing with the call", y="segments") +
  theme_minimal(9) + theme(legend.position="none")
suppressWarnings(ggsave("DR/fig/DR14_1_vote_purity.pdf", p1, width=9, height=7))

hr("2. HOW FAR IS EACH PROPAGATED GENE FROM A MEASURED ONE?")
if (!file.exists("DR/out/DR02_propagated_blocks.csv")) {
  cat("  DR02_propagated.csv not found -- skipping\n")
} else {
  PR <- read_csv("DR/out/DR02_propagated_blocks.csv", show_col_types=FALSE)
  cat(sprintf("  propagated genes: %d\n\n", nrow(PR)))
  votes_by_reg <- GA %>% transmute(genome, chr, region, vmb=mb)
  DIST <- PR %>% mutate(mb=mid/1e6) %>%
    inner_join(votes_by_reg, by=c("genome","chr","region"),
               relationship="many-to-many") %>%
    group_by(gene, genome, chr, region, label) %>%
    summarise(mb=mb[1], d_nearest=min(abs(mb - vmb)), .groups="drop")
  print(as.data.frame(DIST %>% group_by(genome) %>%
    summarise(genes=n(),
              median_kb=round(1000*median(d_nearest),1),
              p75_kb=round(1000*quantile(d_nearest,.75),1),
              p95_Mb=round(quantile(d_nearest,.95),2),
              max_Mb=round(max(d_nearest),2),
              within_100kb=round(mean(d_nearest < 0.1),3),
              beyond_5Mb=round(mean(d_nearest > 5),3), .groups="drop")),
    row.names=FALSE)
  cat("\n  within_100kb is the well-supported fraction. beyond_5Mb is\n")
  cat("  extrapolation and should be drawn faded, or excluded.\n")
  write_csv(DIST, "DR/out/DR14_gene_confidence.csv")
  p2 <- ggplot(DIST, aes(d_nearest)) +
    geom_histogram(bins=60, fill="#378ADD", colour=NA) +
    geom_vline(xintercept=c(0.1, 5), linetype="dashed", colour="firebrick") +
    facet_wrap(~genome, scales="free_y", ncol=2) + scale_x_log10() +
    labs(title="DR14.2 - distance from a propagated gene to the nearest measured one",
         subtitle="dashed at 100 kb and 5 Mb | log scale",
         x="distance to nearest voting gene (Mb)", y="genes") +
    theme_minimal(9)
  suppressWarnings(ggsave("DR/fig/DR14_2_distance.pdf", p2, width=9, height=7))
}

hr("3. DO SEGMENT SPLITS FALL ON SYNTENY BOUNDARIES?")
cat("  The ancestral REGION comes from GENESPACE. The splits WITHIN a region\n")
cat("  came from a changepoint scan on delta: 179 tracks -> 207 segments.\n")
cat("  This asks whether those 28 splits coincide with block edges.\n\n")
nseg <- SS %>% count(genome, region, chr, name="segs")
cat(sprintf("  tracks with >1 segment: %d of %d\n",
            sum(nseg$segs > 1), nrow(nseg)))
if (sum(nseg$segs > 1) == 0) {
  cat("  no splits to check\n")
} else {
  blk <- read_csv(file.path(GSD,"syntenicBlock_coordinates.csv"), show_col_types=FALSE)
  DROS <- paste0("Drosera_", c("regia","binata","paradoxa","scorpioides","capensis"))
  edges <- bind_rows(
    blk %>% filter(genome1 %in% DROS) %>%
      transmute(genome=genome1, chr=chr1,
                e1=pmin(startBp1,endBp1)/1e6, e2=pmax(startBp1,endBp1)/1e6),
    blk %>% filter(genome2 %in% DROS) %>%
      transmute(genome=genome2, chr=chr2,
                e1=pmin(startBp2,endBp2)/1e6, e2=pmax(startBp2,endBp2)/1e6)) %>%
    distinct() %>% pivot_longer(c(e1,e2), values_to="edge") %>% select(-name)
  cat(sprintf("  GENESPACE block edges on Drosera chromosomes: %d\n\n", nrow(edges)))
  brk <- SS %>% group_by(genome, region, chr) %>% filter(n() > 1) %>%
    arrange(mb_lo, .by_group=TRUE) %>%
    summarise(breaks = list(head(mb_hi, -1)), .groups="drop") %>%
    tidyr::unnest(breaks)
  brk$d_edge <- mapply(function(g, c, b) {
    e <- edges$edge[edges$genome==g & edges$chr==c]
    if (!length(e)) NA_real_ else min(abs(e - b))
  }, brk$genome, brk$chr, brk$breaks)
  print(as.data.frame(brk %>% transmute(genome, region, chr,
    break_mb=round(breaks,2), dist_to_block_edge_Mb=round(d_edge,3))),
    row.names=FALSE)
  set.seed(1)
  nullb <- replicate(999, {
    d <- vapply(seq_len(nrow(brk)), function(i) {
      e <- edges$edge[edges$genome==brk$genome[i] & edges$chr==brk$chr[i]]
      if (!length(e)) return(NA_real_)
      rb <- runif(1, min(e), max(e)); min(abs(e - rb))
    }, numeric(1))
    median(d, na.rm=TRUE)
  })
  obs <- median(brk$d_edge, na.rm=TRUE)
  cat(sprintf("\n  median distance to nearest block edge: %.3f Mb\n", obs))
  cat(sprintf("  random-position null: %.3f Mb [%.3f, %.3f]\n",
              median(nullb), quantile(nullb,.05), quantile(nullb,.95)))
  cat(sprintf("  p = %.3f\n", mean(nullb <= obs)))
  cat("\n  observed CLOSER than the null -> splits track synteny, independent\n")
  cat("  corroboration. Same as the null -> they follow delta noise.\n")
}

hr("VERDICT")
cat("  1. purity   -- how many segments rest on a clean vote vs a marginal one\n")
cat("  2. distance -- how much of the 82-94%% coverage is measured vs interpolated\n")
cat("  3. boundaries -- whether the 28 delta-drawn splits have external support\n")
