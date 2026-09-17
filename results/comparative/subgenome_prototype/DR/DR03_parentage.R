#!/usr/bin/env Rscript
# ============================================================================
# DR03 — DID ALL DROSERA INHERIT THE SAME TWO PROGENITORS?
#
# THE MEASUREMENT
#   Two species, one locus, three depths using DR02's segment labels:
#     A<->A   how different are their A copies
#     B<->B   how different are their B copies
#     A<->B   how different is A from B  <- the reference, the A/B split
#
#   Divide the first two by the third. Both now sit on a scale where
#     0 = "diverged only at speciation"
#     1 = "as different as A is from B"
#
#   INTERNALLY CALIBRATED, so no rate correction: numerator and denominator
#   are the same two species at the same locus, so lineage rates and the
#   gene's own rate divide out.
#
# THE THREE OUTCOMES
#   both low          -> shared A and shared B, one hybridisation
#   A high, B low     -> different A progenitor, shared B  (the mosaic case)
#   both high         -> separate hybridisation entirely
#
#   WHY BOTH AXES: the ratio A<->A / B<->B alone would be BLIND to a fully
#   separate hybridisation, because if both came from elsewhere that ratio
#   returns to 1 and looks shared. The A<->B denominator supplies the scale.
#
#   Expected under shared parentage: ~0.35 for close species (speciation
#   T~0.20 over A/B split T~0.57); ~0.70 for regia, which is more distant.
#   Still well below 1 in both cases.
#
# WHAT THIS CANNOT DO
#   It gives DEPTHS, not topology. Whether regia is sister to Drosera or to
#   Dionaea is DR05's tree, not this.
#
# HANDLING THE TWO AWKWARD SPECIES
#   capensis is 12-ploid, probably AAAABB: FOUR A-derived subgenomes against
#   everyone else's two. Taking a minimum over more copies biases it low, so
#   capensis would look artificially close to everything. Its recent-WGD
#   duplicates are collapsed first at dS < 0.25 (script 15's threshold), one
#   representative per cluster, which restores comparability.
#
#   paradoxa carries converted copies (DR01: d_sis p05 = 0.000). A converted
#   copy is artificially close to its partner and would shrink A<->A,
#   mimicking shared parentage. Flagged loci excluded.
#
# EXCLUSIONS
#   - AB03 exchange tracts: Dionaea's labels invert there, so every Drosera
#     label inherited from them inverts too.
#   - loci where a species has no confidently labelled copy of the subgenome
#     being compared.
#
# IN   DR/out/DR02_gene_delta.csv, DR/out/DR02_segments.csv,
#      DR/out/pairwise_ks.csv, DR/out/DR01_votes.csv,
#      AB/out/AB03_flagged_genes.csv
# OUT  DR/out/DR03_locus_ratios.csv   per locus per pair, all three depths
#      DR/out/DR03_summary.csv        the ten points with bootstrap CIs
# FIG  DR/fig/DR03_1_scatter.pdf      THE RESULT - ten pairs, one plot
#      DR/fig/DR03_2_distributions.pdf  the ratios behind each point
#      DR/fig/DR03_3_depths.pdf       raw depths, saturation check
#      DR/fig/DR03_4_control.pdf      label-shuffle null
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))
set.seed(1); NBOOT <- 2000L; COLLAPSE <- 0.25
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

hr("0. LABELS from DR02 segments")
G  <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE)
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE)
lab <- SS %>% filter(label %in% c("A","B")) %>%
       select(genome, region, chr, mb_lo, mb_hi, label)
GL <- G %>% mutate(mb=mid/1e6) %>%
  inner_join(lab, by=c("genome","region","chr"), relationship="many-to-many") %>%
  filter(mb>=mb_lo, mb<=mb_hi) %>% group_by(tip) %>% slice(1) %>% ungroup()
cat(sprintf("  genes inheriting a segment label: %d of %d (%.0f%%)\n",
            nrow(GL), nrow(G), 100*nrow(GL)/nrow(G)))

hr("1. EXCLUSIONS")
n0 <- nrow(GL)
ff <- "AB/out/AB03_flagged_genes.csv"
if (file.exists(ff)) {
  fl <- read_csv(ff, show_col_types=FALSE) %>% filter(in_flip_region)
  CAL <- tibble(pair=paste0("chr",1:8),
                region=c("chr2_dom","chr3_dom","chr4_dom","chr5_dom",
                         "chr6_dom","chr7_dom","chr9_dom","chr10_dom"))
  bad <- fl %>% left_join(CAL, by="pair") %>% distinct(nep_gene, region)
  GL <- GL %>% anti_join(bad, by=c("locus"="nep_gene"))
  cat(sprintf("  AB03 exchange tracts: dropped %d genes\n", n0-nrow(GL)))
} else cat("  AB03 flags not found -- exchange tracts NOT excluded\n")
n1 <- nrow(GL)
if (file.exists("DR/out/DR01_votes.csv")) {
  conv <- read_csv("DR/out/DR01_votes.csv", show_col_types=FALSE) %>%
          filter(genome=="Drosera_paradoxa", d_sis < 0.05) %>% distinct(locus)
  GL <- GL %>% filter(!(genome=="Drosera_paradoxa" & locus %in% conv$locus))
  cat(sprintf("  paradoxa conversion loci (d_sis<0.05): dropped %d genes\n",
              n1-nrow(GL)))
}

hr("2. CAPENSIS — collapse recent-WGD duplicates")
cat("  12-ploid, probably AAAABB: 4 A-derived subgenomes vs everyone else's 2.\n")
cat("  A minimum over more copies is biased low, so capensis would look\n")
cat("  artificially close to everything. Collapse at dS < 0.25 first.\n\n")
K <- read_csv("DR/out/pairwise_ks.csv", show_col_types=FALSE) %>%
     filter(!is.na(dS), dS>=0, dS<5, codons>=100)
DD <- bind_rows(K %>% transmute(locus=anchor, a=seq1, b=seq2, d=dS),
                K %>% transmute(locus=anchor, a=seq2, b=seq1, d=dS))
dmap <- setNames(DD$d, paste(DD$locus, DD$a, DD$b))
gdv <- function(l,x,y) unname(dmap[paste(l,x,y)])
capg <- GL %>% filter(genome=="Drosera_capensis")
cat(sprintf("  capensis copies before: %d\n", nrow(capg)))
keep <- capg %>% group_by(locus, label) %>% group_modify(function(d, key) {
  if (nrow(d) <= 1) return(d)
  n <- nrow(d); cl <- seq_len(n)
  for (i in 1:(n-1)) for (j in (i+1):n) {
    dd <- gdv(key$locus, d$tip[i], d$tip[j])
    if (!is.na(dd) && dd < COLLAPSE) cl[cl==cl[j]] <- cl[i]
  }
  d[!duplicated(cl), ]
}) %>% ungroup()
cat(sprintf("  capensis copies after collapse: %d (%.0f%% kept)\n",
            nrow(keep), 100*nrow(keep)/nrow(capg)))
GL <- bind_rows(GL %>% filter(genome!="Drosera_capensis"), keep)
cat("\n  copies per locus per species (after all handling):\n")
print(as.data.frame(GL %>% count(genome, locus, label) %>%
  group_by(genome, label) %>%
  summarise(median_copies=median(n), max_copies=max(n), .groups="drop") %>%
  pivot_wider(names_from=label, values_from=c(median_copies, max_copies))),
  row.names=FALSE)

hr("3. THE THREE DEPTHS, per locus per species pair")
cat("  A<->A and B<->B: MINIMUM across available copies. All A copies descend\n")
cat("  from the same ancestral A, so under shared parentage the minimum sits\n")
cat("  at speciation depth whichever duplicate pairs with which.\n")
cat("  A<->B: mean of the two directions, both measuring the A/B split.\n")
cat("  Ratios are computed PER LOCUS then medianed -- median of ratios, NOT\n")
cat("  ratio of medians (the script 14 bug).\n\n")
SP <- sort(unique(GL$genome))
mind <- function(l, t1, t2) {
  v <- unlist(lapply(t1, function(x) vapply(t2, function(y) gdv(l,x,y), numeric(1))))
  v <- v[!is.na(v)]; if (!length(v)) NA_real_ else min(v)
}
LR <- bind_rows(lapply(seq_along(SP), function(i)
  bind_rows(lapply(SP[-(1:i)], function(s2) {
    s1 <- SP[i]
    d1 <- GL %>% filter(genome==s1); d2 <- GL %>% filter(genome==s2)
    loci <- intersect(d1$locus, d2$locus)
    bind_rows(lapply(loci, function(l) {
      a1 <- d1$tip[d1$locus==l & d1$label=="A"]; b1 <- d1$tip[d1$locus==l & d1$label=="B"]
      a2 <- d2$tip[d2$locus==l & d2$label=="A"]; b2 <- d2$tip[d2$locus==l & d2$label=="B"]
      if (!length(a1)||!length(a2)||!length(b1)||!length(b2)) return(NULL)
      AA <- mind(l,a1,a2); BB <- mind(l,b1,b2)
      AB <- mean(c(mind(l,a1,b2), mind(l,b1,a2)), na.rm=TRUE)
      if (any(is.na(c(AA,BB,AB))) || AB<=0) return(NULL)
      tibble(pair=paste(sub("Drosera_","",s1), sub("Drosera_","",s2), sep="-"),
             locus=l, AA=AA, BB=BB, AB=AB, rAA=AA/AB, rBB=BB/AB)
    }))
  }))))
cat(sprintf("  loci with all three depths: %d across %d pairs\n",
            nrow(LR), n_distinct(LR$pair)))
print(as.data.frame(LR %>% count(pair, name="loci") %>% arrange(desc(loci))),
      row.names=FALSE)
write_csv(LR, "DR/out/DR03_locus_ratios.csv")

hr("4. RAW DEPTHS — saturation check before ratios")
print(as.data.frame(LR %>% group_by(pair) %>%
  summarise(n=n(), med_AA=round(median(AA),3), med_BB=round(median(BB),3),
            med_AB=round(median(AB),3), frac_AB_over_2=round(mean(AB>2),3),
            .groups="drop") %>% arrange(med_AB)), row.names=FALSE)
cat("\n  A<->B is the A/B split for that pair; Dionaea gave 0.588 (DR02).\n")
cat("  frac_AB_over_2 flags saturation: dS above ~2 is unreliable.\n")
p3 <- ggplot(LR %>% pivot_longer(c(AA,BB,AB), names_to="cmp", values_to="dS"),
             aes(dS, fill=cmp)) +
  geom_density(alpha=0.45, colour=NA) + facet_wrap(~pair, scales="free_y", ncol=2) +
  xlim(0,2.5) +
  scale_fill_manual(values=c(AA="#1D9E75", BB="#D85A30", AB="grey55")) +
  labs(title="DR03.3 - raw depths per species pair",
       subtitle="AA and BB should sit LEFT of AB if both subgenomes are shared",
       x="dS", y="density", fill=NULL) +
  theme_minimal(8) + theme(legend.position="top")
suppressWarnings(ggsave("DR/fig/DR03_3_depths.pdf", p3, width=9, height=11))

hr("5. THE RATIOS — distributions before medians")
cat("  A single median would hide a bimodal ratio distribution, which is what\n")
cat("  LOCUS-LEVEL mosaicism looks like: some loci shared, others not.\n\n")
print(as.data.frame(LR %>% group_by(pair) %>%
  summarise(n=n(),
            rAA_q1=round(quantile(rAA,.25),3), rAA_med=round(median(rAA),3),
            rAA_q3=round(quantile(rAA,.75),3),
            rBB_q1=round(quantile(rBB,.25),3), rBB_med=round(median(rBB),3),
            rBB_q3=round(quantile(rBB,.75),3), .groups="drop") %>%
  arrange(rAA_med)), row.names=FALSE)
p2 <- ggplot(LR %>% pivot_longer(c(rAA,rBB), names_to="which", values_to="r"),
             aes(r, fill=which)) +
  geom_density(alpha=0.5, colour=NA) +
  geom_vline(xintercept=1, colour="grey30", linetype="dashed") +
  facet_wrap(~pair, scales="free_y", ncol=2) + xlim(0,2) +
  scale_fill_manual(values=c(rAA="#1D9E75", rBB="#D85A30"),
                    labels=c("A<->A / A<->B","B<->B / A<->B")) +
  labs(title="DR03.2 - scaled divergence distributions",
       subtitle="dashed 1.0 = as divergent as A is from B | bimodal = locus-level mosaicism",
       x="ratio", y="density", fill=NULL) +
  theme_minimal(8) + theme(legend.position="top")
suppressWarnings(ggsave("DR/fig/DR03_2_distributions.pdf", p2, width=9, height=11))

hr("6. THE TEN POINTS, with bootstrap intervals")
bmed <- function(x) { b <- replicate(NBOOT, median(sample(x, replace=TRUE)))
                      quantile(b, c(.025,.975)) }
SUM <- LR %>% group_by(pair) %>%
  group_modify(function(d, key) {
    ca <- bmed(d$rAA); cb <- bmed(d$rBB)
    tibble(n=nrow(d), rAA=median(d$rAA), rAA_lo=ca[1], rAA_hi=ca[2],
           rBB=median(d$rBB), rBB_lo=cb[1], rBB_hi=cb[2],
           w=suppressWarnings(wilcox.test(d$rAA, d$rBB, paired=TRUE)$p.value))
  }) %>% ungroup() %>%
  mutate(has_regia=grepl("regia", pair),
         diff=round(rAA-rBB,3), p_adj=p.adjust(w,"BH"))
print(as.data.frame(SUM %>% transmute(pair, n,
  rAA=round(rAA,3), rAA_CI=sprintf("[%.2f,%.2f]", rAA_lo, rAA_hi),
  rBB=round(rBB,3), rBB_CI=sprintf("[%.2f,%.2f]", rBB_lo, rBB_hi),
  diff, p_adj=signif(p_adj,3)) %>% arrange(desc(diff))), row.names=FALSE)
cat("\n  diff = rAA - rBB. Positive and significant = A more divergent than B\n")
cat("  for that pair, i.e. a different A progenitor with a shared B.\n")
write_csv(SUM, "DR/out/DR03_summary.csv")
p1 <- ggplot(SUM, aes(rAA, rBB)) +
  geom_abline(slope=1, intercept=0, colour="grey70", linetype="dotted") +
  geom_hline(yintercept=1, colour="grey70", linetype="dashed") +
  geom_vline(xintercept=1, colour="grey70", linetype="dashed") +
  geom_linerange(aes(xmin=rAA_lo, xmax=rAA_hi), colour="grey60") +
  geom_linerange(aes(ymin=rBB_lo, ymax=rBB_hi), colour="grey60") +
  geom_point(aes(colour=has_regia), size=3) +
  ggrepel::geom_text_repel(aes(label=pair), size=2.6, max.overlaps=20) +
  scale_colour_manual(values=c("FALSE"="#378ADD","TRUE"="#D85A30"),
                      labels=c("no regia","involves regia")) +
  coord_equal(xlim=c(0,1.3), ylim=c(0,1.3)) +
  labs(title="DR03.1 - shared parentage per subgenome",
       subtitle="bottom-left = both shared | bottom-right = different A, shared B | top-right = separate hybridisation",
       x="A<->A / A<->B", y="B<->B / A<->B", colour=NULL) +
  theme_minimal(11) + theme(legend.position="top")
suppressWarnings(ggsave("DR/fig/DR03_1_scatter.pdf", p1, width=8, height=8))

hr("7. CONTROL — shuffle the labels")
cat("  If A/B labels carry no information, A<->A and A<->B measure the same\n")
cat("  thing and BOTH ratios go to 1. Observed ratios well below the shuffled\n")
cat("  null mean the labels are doing real work.\n\n")
SHUF <- bind_rows(lapply(unique(LR$pair), function(pp) {
  d <- LR %>% filter(pair==pp)
  tibble(pair=pp, r=c(d$AA, d$BB)/d$AB[c(seq_len(nrow(d)), seq_len(nrow(d)))],
         src="observed")
}))
nullr <- bind_rows(lapply(unique(GL$genome), function(g) NULL))
cat("  (label-shuffle null: recompute with labels permuted within each locus)\n")
GS <- GL %>% group_by(locus, genome) %>% mutate(label=sample(label)) %>% ungroup()
LRS <- bind_rows(lapply(seq_along(SP), function(i)
  bind_rows(lapply(SP[-(1:i)], function(s2) {
    s1 <- SP[i]
    d1 <- GS %>% filter(genome==s1); d2 <- GS %>% filter(genome==s2)
    loci <- intersect(d1$locus, d2$locus)
    bind_rows(lapply(loci, function(l) {
      a1 <- d1$tip[d1$locus==l & d1$label=="A"]; b1 <- d1$tip[d1$locus==l & d1$label=="B"]
      a2 <- d2$tip[d2$locus==l & d2$label=="A"]; b2 <- d2$tip[d2$locus==l & d2$label=="B"]
      if (!length(a1)||!length(a2)||!length(b1)||!length(b2)) return(NULL)
      AA <- mind(l,a1,a2); AB <- mean(c(mind(l,a1,b2), mind(l,b1,a2)), na.rm=TRUE)
      if (any(is.na(c(AA,AB))) || AB<=0) return(NULL)
      tibble(pair=paste(sub("Drosera_","",s1), sub("Drosera_","",s2), sep="-"),
             rAA=AA/AB)
    }))
  }))))
cmp <- LR %>% group_by(pair) %>% summarise(obs=median(rAA), .groups="drop") %>%
  left_join(LRS %>% group_by(pair) %>%
              summarise(shuffled=median(rAA), .groups="drop"), by="pair")
cmp$gap <- round(cmp$shuffled - cmp$obs, 3)
print(as.data.frame(cmp %>% mutate(obs=round(obs,3), shuffled=round(shuffled,3)) %>%
  arrange(desc(gap))), row.names=FALSE)
cat("\n  gap > 0 means real labels give a LOWER A<->A than shuffled ones,\n")
cat("  i.e. the labels are identifying genuinely closer copies.\n")
p4 <- ggplot(bind_rows(LR %>% transmute(pair, rAA, src="observed"),
                       LRS %>% transmute(pair, rAA, src="labels shuffled")),
             aes(rAA, fill=src)) +
  geom_density(alpha=0.5, colour=NA) +
  facet_wrap(~pair, scales="free_y", ncol=2) + xlim(0,2) +
  scale_fill_manual(values=c("observed"="#378ADD","labels shuffled"="grey55")) +
  labs(title="DR03.4 - control: does the A/B labelling do any work?",
       subtitle="observed shifted LEFT of shuffled = labels identify genuinely closer copies",
       x="A<->A / A<->B", y="density", fill=NULL) +
  theme_minimal(8) + theme(legend.position="top")
suppressWarnings(ggsave("DR/fig/DR03_4_control.pdf", p4, width=9, height=11))

hr("READ IN THIS ORDER")
cat("  1. sec 7  -- did the labels beat a shuffle? if not, stop.\n")
cat("  2. sec 4  -- are the raw depths sane? A<->B near 0.6-1.2, not saturated.\n")
cat("  3. FIG .2 -- are the ratio distributions unimodal? bimodal = mosaic loci.\n")
cat("  4. FIG .1 -- THE RESULT. where do the ten pairs land?\n")
cat("  5. sec 6  -- is rAA - rBB significant for any regia pair?\n")
