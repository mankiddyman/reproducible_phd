#!/usr/bin/env Rscript
# ============================================================================
# DR01 — IS ANY DROSERA SPECIES ABC RATHER THAN AAB?   [v2]
#
# THE QUESTION
#   A hexaploid has three subgenomes. Two ancestral lineages (AAB) or three
#   (ABC)? Under AAB two copies share a more recent ancestor; under ABC all
#   three are equidistant.
#
# THE TEST
#   Four-point condition (Buneman 1971) on the three copies + Nepenthes.
#   Three ways to split four tips into two pairs; smallest distance sum wins
#   and names the sister pair. RATE-FREE: each tip appears exactly once in
#   each sum, so lineage rate cancels from the ranking.
#   Per region, AAB predicts the SAME pair wins repeatedly. Chi-square vs
#   uniform, 2 df -- not "max count vs binomial 1/3", which inflates.
#
# ----------------------------------------------------------------------------
# WHAT v1 GOT WRONG, AND WHAT CHANGED
# ----------------------------------------------------------------------------
#  1. THE WINNING PAIR IS THE CLOSEST PAIR, NOT NECESSARILY THE ANCESTRAL ONE.
#     If a species carries A, A' (recent duplicate or converted copy) and B,
#     the four-point picks A/A' and reports "AAB structure" for the wrong
#     reason. v1 read paradoxa's ratio of 0.103 as AAB when d_sis there is
#     ~0.15 against an A/B depth of 1.49 -- a tenth. That is a recent
#     duplication or gene conversion, not the ancestral A1/A2 split.
#     regia by contrast gave d_sis ~0.22 against A/B 0.344, matching the
#     independently established A1|A2 = 0.22 vs A|B = 0.358.
#     FIX: section 3 plots the d_sis distribution FIRST and the chi-square is
#     reported both pooled and stratified.
#
#  2. THE NULL USED ADDITIVE NOISE WHERE NOISE SCALES WITH DISTANCE.
#     v1 estimated one absolute sigma per species. regia's sigma of 0.2256
#     against a median d_out of 0.344 is 65% relative noise, so simulated
#     triples were noise-dominated, the null ratio collapsed to 0.391 and
#     regia looked WORSE than ABC. FIX: relative noise, multiplicative
#     (lognormal) simulation.
#
#  3. HALF THE LOCI SILENTLY VANISHED. 799 of 1524. The dS<3 ceiling kills
#     Nepenthes comparisons, which run ~1.2 median and higher for fast genes.
#     FIX: ceiling raised to 5, and section 1 reports exactly what each
#     filter drops, by pair type.
#
#  4. CHROMOSOME NAMES DO NOT IDENTIFY SUBGENOMES. "chrX|chrX" was treated
#     as a vote category, but if A1 and B share a chromosome while A2 sits
#     elsewhere, that label is identical to the case where A1 and A2 share
#     one. FIX: the pair is named by GENESPACE SYNTENIC BLOCK ID against
#     Nepenthes. Two loci in the same block descend from the same segment,
#     which is what "the same pair wins" actually requires.
#
#  5. dplyr SEQUENTIAL-EVALUATION BUG, 7th recorded instance in this project.
#     summarise(same_chr=sum(same_chr), pct=mean(same_chr)) collapsed the
#     column before mean() saw it, printing 900%. FIX: every derived quantity
#     computed outside the pipeline.
#
#  6. Loci with several Nepenthes tips caused a many-to-many join. FIX: one
#     outgroup tip per locus.
#
# IN   DR/locus_meta.tsv, DR/out/pairwise_ks.csv,
#      ../genespace/results/syntenicBlock_coordinates.csv
# OUT  DR/out/DR01_votes.csv, DR01_regions.csv, DR01_power.csv
# FIG  DR/fig/DR01_0_dsis.pdf     READ FIRST - recent vs ancestral sisters
#      DR/fig/DR01_1_votes.pdf    wins per region
#      DR/fig/DR01_2_gap.pdf      is there a clear winner
#      DR/fig/DR01_3_depth.pdf    depth ratio vs star null
#      DR/fig/DR01_4_power.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))
set.seed(1)
DSMAX <- 5; GAPTHR <- 0.02
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")
GSD <- file.path(dirname(getwd()), "genespace", "results")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

hr("0. ASSEMBLE")
LOC <- read_tsv("DR/locus_meta.tsv", show_col_types=FALSE) %>%
       group_by(tip) %>% slice(1) %>% ungroup()
KR  <- read_csv("DR/out/pairwise_ks.csv", show_col_types=FALSE)
cat(sprintf("  raw dS rows: %d\n", nrow(KR)))
K <- KR %>% filter(!is.na(dS), dS>=0, dS<DSMAX, codons>=100)
cat(sprintf("  after 0<=dS<%d & codons>=100: %d (%.1f%% kept)\n",
            DSMAX, nrow(K), 100*nrow(K)/nrow(KR)))
cat("\n  what the dS ceiling costs, by comparison type:\n")
KR2 <- KR %>% filter(!is.na(dS), dS>=0) %>%
  mutate(type=ifelse(sp1=="Nepenthes_gracilis" | sp2=="Nepenthes_gracilis",
                     "vs Nepenthes", "ingroup"))
drop <- KR2 %>% group_by(type) %>%
  summarise(pairs=n(), median_dS=round(median(dS),3),
            over3=sum(dS>=3), over5=sum(dS>=DSMAX), .groups="drop")
drop$pct_lost_at_3 <- round(100*drop$over3/drop$pairs, 1)
drop$pct_lost_at_5 <- round(100*drop$over5/drop$pairs, 1)
print(as.data.frame(drop), row.names=FALSE)
cat(sprintf("\n  pairs with dS EXACTLY 0 (identical sequences): %d\n",
            sum(KR$dS==0, na.rm=TRUE)))
cat("  these are very recent duplicates or complete conversions -- KEPT,\n")
cat("  because excluding them would hide the class FIG 0 exists to show.\n")

nept <- LOC %>% filter(genome=="Nepenthes_gracilis") %>%
        group_by(locus) %>% slice(1) %>% ungroup() %>%
        transmute(locus, nep_tip=tip)
sel <- LOC %>% filter(genome %in% DROS) %>%
  group_by(locus, genome) %>% filter(n()==3) %>% ungroup() %>%
  inner_join(nept, by="locus") %>% mutate(region=sub("-.*$","",locus))
cat("\n  candidate (locus, species) with exactly 3 copies + outgroup:\n")
print(as.data.frame(sel %>% distinct(locus, genome) %>% count(genome, name="loci")),
      row.names=FALSE)

hr("1. SYNTENIC BLOCK IDs — the pair identifier")
cat("  Chromosome names do not identify subgenomes: two copies can share a\n")
cat("  chromosome after a fusion without being homeologs of each other.\n")
cat("  GENESPACE block ID against Nepenthes is the stable segment identity.\n\n")
blk <- read_csv(file.path(GSD,"syntenicBlock_coordinates.csv"), show_col_types=FALSE)
nb <- bind_rows(
  blk %>% filter(genome2=="Nepenthes_gracilis", genome1 %in% DROS) %>%
    transmute(sp=genome1, blkID, chr=chr1,
              lo=pmin(startBp1,endBp1), hi=pmax(startBp1,endBp1), nh=nHits1),
  blk %>% filter(genome1=="Nepenthes_gracilis", genome2 %in% DROS) %>%
    transmute(sp=genome2, blkID, chr=chr2,
              lo=pmin(startBp2,endBp2), hi=pmax(startBp2,endBp2), nh=nHits2)) %>%
  filter(!is.na(lo), !is.na(hi)) %>% distinct()
cat(sprintf("  Drosera-vs-Nepenthes blocks: %d\n", nrow(nb)))
bedpos <- read.table(file.path(GSD,"combBed.txt"), header=TRUE, sep="\t",
                     quote="", comment.char="", stringsAsFactors=FALSE) %>%
  filter(genome %in% DROS) %>% transmute(genome, gene=id, mid=(start+end)/2)
sel <- sel %>% left_join(bedpos, by=c("genome","gene"))
find_blk <- function(sp_, chr_, mid_) {
  if (is.na(mid_)) return(NA_character_)
  c <- nb[nb$sp==sp_ & nb$chr==chr_ & nb$lo<=mid_ & nb$hi>=mid_, ]
  if (!nrow(c)) return(NA_character_)
  as.character(c$blkID[which.max(c$nh)])
}
sel$blk <- mapply(find_blk, sel$genome, sel$chr, sel$mid)
nblk <- sum(!is.na(sel$blk))
cat(sprintf("  copies assigned to a block: %d of %d (%.1f%%)\n",
            nblk, nrow(sel), 100*nblk/nrow(sel)))
cat("  copies with no block fall back to their chromosome name.\n")
sel <- sel %>% mutate(unit = ifelse(is.na(blk), paste0("chr:", chr), paste0("b", blk)))

hr("2. FOUR-POINT, per locus")
DD <- bind_rows(K %>% transmute(locus=anchor, a=seq1, b=seq2, d=dS),
                K %>% transmute(locus=anchor, a=seq2, b=seq1, d=dS))
dmap <- setNames(DD$d, paste(DD$locus, DD$a, DD$b))
gd <- function(l,x,y) unname(dmap[paste(l,x,y)])
keys <- sel %>% distinct(locus, genome)
miss_reason <- c(ingroup=0L, outgroup=0L)
res <- bind_rows(lapply(seq_len(nrow(keys)), function(i) {
  l <- keys$locus[i]; g <- keys$genome[i]
  d <- sel[sel$locus==l & sel$genome==g, ]
  N <- d$nep_tip[1]; t <- d$tip; u <- d$unit
  m <- matrix(NA_real_, 3, 3)
  for (x in 1:2) for (y in (x+1):3) { m[x,y] <- m[y,x] <- gd(l, t[x], t[y]) }
  dn <- vapply(1:3, function(x) gd(l, t[x], N), numeric(1))
  if (any(is.na(m[upper.tri(m)]))) { miss_reason["ingroup"] <<- miss_reason["ingroup"]+1L; return(NULL) }
  if (any(is.na(dn)))              { miss_reason["outgroup"] <<- miss_reason["outgroup"]+1L; return(NULL) }
  S <- c(m[1,2]+dn[3], m[1,3]+dn[2], m[2,3]+dn[1])
  pr <- list(c(1,2), c(1,3), c(2,3))[[which.min(S)]]
  og <- setdiff(1:3, pr); o <- sort(S)
  tibble(locus=l, genome=g, region=sub("-.*$","",l),
         pair=paste(sort(u[pr]), collapse="|"),
         same_unit=u[pr[1]]==u[pr[2]],
         d_sis=m[pr[1],pr[2]], d_out1=m[pr[1],og], d_out2=m[pr[2],og],
         gap=(o[2]-o[1])/mean(S), excess=(o[3]-o[2])/mean(S))
}))
cat(sprintf("  complete: %d of %d | dropped: %d ingroup, %d outgroup distance missing\n",
            nrow(res), nrow(keys), miss_reason["ingroup"], miss_reason["outgroup"]))
print(as.data.frame(res %>% count(genome, name="loci")), row.names=FALSE)
su <- res %>% count(genome, same_unit) %>%
      pivot_wider(names_from=same_unit, values_from=n, values_fill=0)
cat("\n  wins where both sisters are in the SAME syntenic block:\n")
print(as.data.frame(su), row.names=FALSE)
cat("  (same block = a within-segment duplicate, not a homeolog pair)\n")

hr("3. d_sis DISTRIBUTION — READ THIS FIRST")
cat("  The four-point picks the CLOSEST pair. If a species carries a recent\n")
cat("  duplicate or a converted copy, that pair wins and the test reports\n")
cat("  'AAB' for the wrong reason. Recent duplicates and ancestral A1/A2 sit\n")
cat("  at different depths, so a BIMODAL d_sis distribution reveals the mix.\n\n")
ds <- res %>% group_by(genome) %>%
  summarise(loci=n(), p05=round(quantile(d_sis,.05),3),
            p25=round(quantile(d_sis,.25),3), median=round(median(d_sis),3),
            p75=round(quantile(d_sis,.75),3),
            med_dout=round(median((d_out1+d_out2)/2),3), .groups="drop")
ds$ratio_of_medians <- round(ds$median/ds$med_dout, 3)
print(as.data.frame(ds), row.names=FALSE)
cat("\n  For reference (script 18/23, independently established):\n")
cat("    regia    A1|A2 = 0.22, A|B = 0.358  -> expect d_sis ~0.22\n")
cat("    capensis recent WGD at dS ~0.106    -> a low mode is EXPECTED there\n")
cat("    core three A|B ~1.30                -> ancestral A1|A2 should be near it\n")
p0 <- ggplot(res, aes(d_sis)) +
  geom_histogram(bins=60, fill="#378ADD", colour="white", linewidth=0.2) +
  geom_vline(xintercept=0.106, colour="firebrick", linetype="dashed") +
  geom_vline(xintercept=0.22,  colour="darkgreen", linetype="dashed") +
  facet_wrap(~genome, scales="free", ncol=2) +
  labs(title="DR01.0 - depth of the winning sister pair",
       subtitle="red 0.106 = capensis recent WGD | green 0.22 = regia ancestral A1/A2 | bimodal = a mixture",
       x="d_sis", y="loci") + theme_minimal(9)
suppressWarnings(ggsave("DR/fig/DR01_0_dsis.pdf", p0, width=9, height=8))

hr("4. ADDITIVITY")
ad <- res %>% group_by(genome) %>%
  summarise(median_gap=round(median(gap),4),
            median_excess=round(median(excess),4),
            f01=round(mean(gap>0.01),3), f02=round(mean(gap>0.02),3),
            f05=round(mean(gap>0.05),3), .groups="drop")
print(as.data.frame(ad), row.names=FALSE)
cat("\n  Buneman: for additive distances the two largest sums are EQUAL, so\n")
cat("  excess ~ 0. Non-zero excess = additivity violated, treat with care.\n")
res <- res %>% mutate(clear = gap > GAPTHR)
p2 <- ggplot(res, aes(gap)) +
  geom_histogram(bins=60, fill="#378ADD", colour="white", linewidth=0.2) +
  geom_vline(xintercept=GAPTHR, colour="firebrick", linetype="dashed") +
  facet_wrap(~genome, scales="free_y", ncol=2) + xlim(0,0.4) +
  labs(title="DR01.2 - separation between best and second-best pairing",
       subtitle="mass at zero = near-ties | dashed = clear-winner cut",
       x="gap = (S_mid - S_min) / mean(S)", y="loci") + theme_minimal(9)
suppressWarnings(ggsave("DR/fig/DR01_2_gap.pdf", p2, width=9, height=8))

hr("5. CHI-SQUARE per region — pooled and stratified")
cat("  Stratified excludes loci whose sisters are shallower than 40% of the\n")
cat("  A/B depth, i.e. likely recent duplicates or conversions rather than\n")
cat("  the ancestral A1/A2 split.\n\n")
res <- res %>% mutate(d_out=(d_out1+d_out2)/2,
                      ratio=ifelse(d_out>0, d_sis/d_out, NA_real_),
                      stratum=ifelse(is.na(ratio), "undefined",
                              ifelse(ratio < 0.40, "recent", "ancestral")))
st <- res %>% count(genome, stratum) %>%
      pivot_wider(names_from=stratum, values_from=n, values_fill=0)
print(as.data.frame(st), row.names=FALSE)
chi_one <- function(v) { n <- sum(v); e <- n/3; sum((v-e)^2/e) }
region_test <- function(d) {
  d %>% group_by(genome, region) %>%
    group_modify(function(x, key) {
      tb <- sort(table(x$pair), decreasing=TRUE)
      v <- as.integer(tb); v <- c(v, rep(0L, max(0, 3-length(v))))[1:3]
      tibble(n=nrow(x), n_units=length(tb), top_pair=names(tb)[1],
             top_n=v[1], chi2=chi_one(v))
    }) %>% ungroup() %>%
    mutate(top_frac=round(top_n/n,3), p=1-pchisq(chi2, df=2)) %>%
    group_by(genome) %>% mutate(p_adj=p.adjust(p,"BH")) %>% ungroup()
}
REG_all <- region_test(res %>% filter(clear)) %>% mutate(set="pooled")
REG_anc <- region_test(res %>% filter(clear, stratum=="ancestral")) %>%
           mutate(set="ancestral only")
REG <- bind_rows(REG_all, REG_anc)
for (s in c("pooled","ancestral only")) {
  cat(sprintf("\n  --- %s ---\n", s))
  print(as.data.frame(REG %>% filter(set==s, n>=15) %>%
    transmute(genome, region, n, n_units, top_frac,
              chi2=round(chi2,2), p_adj=signif(p_adj,3)) %>%
    arrange(genome, desc(n))), row.names=FALSE)
}
cat("\n  n_units > 3 means the region spans more than three segments, so the\n")
cat("  1/3 null understates the true number of categories -- read with care.\n")
write_csv(REG, "DR/out/DR01_regions.csv")
plotd <- res %>% filter(clear) %>% group_by(genome, region) %>%
         mutate(nn=n()) %>% filter(nn>=15) %>% ungroup()
p1 <- ggplot(plotd, aes(region, fill=pair)) +
  geom_bar(position="fill") +
  geom_hline(yintercept=c(1/3,2/3), colour="grey30", linetype="dotted") +
  facet_wrap(~genome, scales="free_x", ncol=2) +
  labs(title="DR01.1 - which syntenic-block pair is sister, per region",
       subtitle="AAB = one colour dominates | ABC = even thirds | dotted = 1/3, 2/3",
       x="ancestral region", y="fraction of loci") +
  theme_minimal(9) + theme(legend.position="none")
suppressWarnings(ggsave("DR/fig/DR01_1_votes.pdf", p1, width=10, height=8))

hr("6. DEPTH RATIO vs a MULTIPLICATIVE star null")
cat("  v1 used ADDITIVE noise. regia's sigma was 65% of its median distance,\n")
cat("  so simulated triples were noise-dominated and the null collapsed.\n")
cat("  Noise in dS scales with dS, so it is estimated and applied RELATIVELY.\n\n")
DR <- res %>% filter(clear, !is.na(ratio))
rel <- DR %>% group_by(genome) %>%
  summarise(rel_sd = sd((d_out1-d_out2)/((d_out1+d_out2)/2), na.rm=TRUE)/sqrt(2),
            .groups="drop")
print(as.data.frame(rel %>% mutate(rel_sd=round(rel_sd,3))), row.names=FALSE)
cat("  rel_sd is relative measurement error on a single dS. Above ~0.4 the\n")
cat("  two sisters plainly differ in rate and the null is conservative.\n\n")
NULLD <- bind_rows(lapply(unique(DR$genome), function(g) {
  d <- DR %>% filter(genome==g); s <- rel$rel_sd[rel$genome==g]
  mu <- rowMeans(cbind(d$d_sis, d$d_out1, d$d_out2))
  x <- matrix(rep(mu,3) * exp(rnorm(3*length(mu), 0, s)), ncol=3)
  tibble(genome=g, ratio=apply(x,1,min)/apply(x,1,function(v) mean(sort(v)[2:3])),
         src="star null (ABC)")
}))
obs <- DR %>% group_by(genome) %>%
  summarise(loci=n(), obs_med=round(median(ratio),3), .groups="drop")
nul <- NULLD %>% group_by(genome) %>%
  summarise(null_med=round(median(ratio),3),
            null_p05=round(quantile(ratio,.05),3), .groups="drop")
cmp <- obs %>% left_join(nul, by="genome")
cmp$shift <- round(cmp$obs_med - cmp$null_med, 3)
cmp$verdict <- ifelse(cmp$shift < -0.05, "sisters shallower than ABC (AAB)",
               ifelse(cmp$shift >  0.05, "NOT shallower", "indistinguishable"))
print(as.data.frame(cmp), row.names=FALSE)
p3 <- ggplot(bind_rows(DR %>% transmute(genome, ratio, src="observed"), NULLD),
             aes(ratio, fill=src)) +
  geom_histogram(aes(y=after_stat(density)), bins=50, alpha=0.55,
                 position="identity", colour=NA) +
  facet_wrap(~genome, scales="free_y", ncol=2) + xlim(0,1.2) +
  scale_fill_manual(values=c("star null (ABC)"="grey55","observed"="#378ADD")) +
  labs(title="DR01.3 - depth ratio: observed vs multiplicative ABC star null",
       subtitle="observed LEFT of null = two copies genuinely closer = AAB",
       x="d_sis / mean(d_out)", y="density", fill=NULL) +
  theme_minimal(9) + theme(legend.position="top")
suppressWarnings(ggsave("DR/fig/DR01_3_depth.pdf", p3, width=9, height=8))
write_csv(res, "DR/out/DR01_votes.csv")

hr("7. POWER")
POW <- bind_rows(lapply(c(15,25,40,60,100,200), function(n)
  bind_rows(lapply(c(0.40,0.50,0.60,0.70), function(f) {
    hit <- mean(replicate(1000, {
      v <- as.integer(table(factor(sample(1:3, n, TRUE,
             prob=c(f,(1-f)/2,(1-f)/2)), levels=1:3)))
      chi_one(v) > qchisq(0.95, 2) }))
    tibble(n=n, true_frac=f, detect=round(hit,2)) }))))
print(as.data.frame(POW %>% pivot_wider(names_from=true_frac, values_from=detect,
                                        names_prefix="frac=")), row.names=FALSE)
write_csv(POW, "DR/out/DR01_power.csv")
p4 <- ggplot(POW, aes(n, detect, colour=factor(true_frac))) +
  geom_hline(yintercept=0.8, linetype="dashed", colour="grey50") +
  geom_line(linewidth=0.8) + geom_point(size=1.6) + ylim(0,1) +
  labs(title="DR01.4 - power of the chi-square test", x="loci per region",
       y="detection rate", colour="true fraction") + theme_minimal(10)
suppressWarnings(ggsave("DR/fig/DR01_4_power.pdf", p4, width=9, height=6))

hr("VERDICT")
va <- REG %>% filter(set=="ancestral only", n>=15)
V <- va %>% group_by(genome) %>%
  summarise(regions=n(), significant=sum(p_adj<0.05),
            med_top_frac=round(median(top_frac),3),
            med_n=median(n), .groups="drop") %>%
  left_join(cmp %>% select(genome, obs_med, null_med, shift), by="genome")
V$verdict <- with(V, ifelse(is.na(shift), "no data",
  ifelse(significant >= 0.6*regions & shift < -0.05, "AAB",
  ifelse(significant == 0 & med_n >= 25, "ABC or no structure",
         "MIXED / underpowered"))))
print(as.data.frame(V), row.names=FALSE)
cat("\n  READ IN THIS ORDER:\n")
cat("   1. FIG .0  -- is d_sis bimodal? that separates recent from ancestral.\n")
cat("   2. sec 5   -- do pooled and ancestral-only disagree? if so, mixture.\n")
cat("   3. FIG .1  -- does one colour dominate per region?\n")
cat("   4. FIG .3  -- is observed left of the null?\n")
cat("   5. sec 7   -- were the non-significant regions powered?\n")
