#!/usr/bin/env Rscript
# ============================================================================
# DR05e — IS regia SISTER TO DIONAEA? THE RATE-FREE TEST
#
# WHY THIS EXISTS
#   Four analyses now place regia with Dionaea:
#     ASTRAL A   (regia_A,Dionaea_A) 0.198 CU -> 0.223 on informative trees
#     ASTRAL B   (regia_B,Dionaea_B) 0.192 CU -> 0.250 on informative trees
#     DR03b      rAA(regia,Dio) 0.649 vs core 0.820, p = 8e-46
#     DR03       rAA ~= rBB for regia -- both subgenomes shift EQUALLY, which
#                is what a shared shift in regia's position predicts
#   The branch GROWING on informative gene trees refutes rate attraction,
#   which predicts the opposite.
#
#   But all four are likelihood- or distance-based, so a rate confound is not
#   fully excluded. Buneman's four-point condition is RATE-FREE: each of the
#   four tips enters each of the three distance sums exactly once, so any
#   lineage's extra branch length adds equally to all three and the RANKING is
#   unchanged. It cannot be fooled by regia being slow.
#
#   And concatenation is measurably wrong at this node. sCF at A_full node 8
#   (the branch separating regia_A) is 32.66 with sDF2 = 41.61 -- at chance,
#   and the ALTERNATIVE resolution has more site support than the tree's own.
#   Bootstrap says 100/100. That is exactly what sCF is for.
#
# THE QUARTET
#   (regia_X, Dionaea_X, CORE_X, Nepenthes)   for X in {A,B}, CORE in the four
#     (regia,Dionaea)|(core,Nep)  -> regia SISTER TO DIONAEA
#     (regia,core)|(Dionaea,Nep)  -> regia INSIDE DROSERA
#     (regia,Nep)|(Dionaea,core)  -> neither; regia outside both
#
# INDEPENDENCE — the thing that would inflate the p-value
#   8 quartets per locus, but they share regia, Dionaea and Nepenthes. Treating
#   them as independent multiplies n by ~8 and gives a spuriously tiny p.
#   So: score every quartet, take the LOCUS majority, test locus majorities
#   against a 1/3 null. Per-quartet counts are reported for description only.
#
# ADDITIVITY
#   Buneman holds exactly for additive distances; real data almost always
#   violates this. Loci where the best and second-best sums are near-tied carry
#   no topological information, so the gap is reported and the test run at
#   several thresholds.
#
# IN   DR/out/DR02_gene_delta.csv, DR/out/DR02_segments.csv,
#      DR/out/pairwise_ks.csv, DR/locus_meta.tsv, fractionation_by_chrpair.csv
# OUT  DR/out/DR05e_quartets.csv, DR05e_loci.csv
# FIG  DR/fig/DR05e_1_votes.pdf, DR05e_2_gap.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))
CORE <- c("binata","paradoxa","scorpioides","capensis")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

hr("1. ASSEMBLE labelled copies")
LOC <- read_tsv("DR/locus_meta.tsv", show_col_types=FALSE) %>%
       group_by(tip) %>% slice(1) %>% ungroup()
K <- read_csv("DR/out/pairwise_ks.csv", show_col_types=FALSE) %>%
     filter(!is.na(dS), dS>=0, dS<5, codons>=100)
DD <- bind_rows(K %>% transmute(locus=anchor, a=seq1, b=seq2, d=dS),
                K %>% transmute(locus=anchor, a=seq2, b=seq1, d=dS))
dmap <- setNames(DD$d, paste(DD$locus, DD$a, DD$b))
gd <- function(l,x,y) ifelse(x==y, 0, unname(dmap[paste(l,x,y)]))
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
G  <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE)
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE)
lab <- SS %>% filter(label %in% c("A","B")) %>%
       select(genome, region, chr, mb_lo, mb_hi, label)
GL <- G %>% mutate(mb=mid/1e6) %>%
  inner_join(lab, by=c("genome","region","chr"), relationship="many-to-many") %>%
  filter(mb>=mb_lo, mb<=mb_hi) %>% group_by(tip) %>% slice(1) %>% ungroup() %>%
  mutate(sp=sub("Drosera_","",genome))
cat(sprintf("  labelled Drosera copies: %d | Dionaea pairs: %d | Nepenthes loci: %d\n",
            nrow(GL), nrow(dio), nrow(nep)))

hr("2. FOUR-POINT, every quartet")
base <- dio %>% inner_join(nep, by="locus")
reg <- GL %>% filter(sp=="regia") %>% select(locus, label, rtip=tip)
cor <- GL %>% filter(sp %in% CORE) %>% select(locus, label, csp=sp, ctip=tip)
Q <- reg %>% inner_join(cor, by=c("locus","label"), relationship="many-to-many") %>%
     inner_join(base, by="locus") %>%
     mutate(Dio = ifelse(label=="A", DioA, DioB))
cat(sprintf("  candidate quartets: %d\n", nrow(Q)))
Q <- Q %>%
  mutate(d_rD = gd(locus, rtip, Dio),  d_cN = gd(locus, ctip, Nep),
         d_rc = gd(locus, rtip, ctip), d_DN = gd(locus, Dio,  Nep),
         d_rN = gd(locus, rtip, Nep),  d_Dc = gd(locus, Dio,  ctip)) %>%
  filter(if_all(starts_with("d_"), ~ !is.na(.))) %>%
  mutate(S1 = d_rD + d_cN, S2 = d_rc + d_DN, S3 = d_rN + d_Dc)
best <- apply(Q[,c("S1","S2","S3")], 1, which.min)
srt  <- t(apply(Q[,c("S1","S2","S3")], 1, sort))
Q <- Q %>% mutate(
  vote = c("regia+Dionaea","regia+core","neither")[best],
  gap  = (srt[,2]-srt[,1]) / rowMeans(cbind(Q$S1,Q$S2,Q$S3)))
cat(sprintf("  complete quartets: %d across %d loci\n", nrow(Q), n_distinct(Q$locus)))
cat("\n  per-quartet votes (DESCRIPTIVE ONLY -- quartets share 3 of 4 tips):\n")
print(as.data.frame(Q %>% count(label, vote) %>%
  pivot_wider(names_from=vote, values_from=n, values_fill=0)), row.names=FALSE)
print(as.data.frame(Q %>% count(csp, label, vote) %>%
  pivot_wider(names_from=vote, values_from=n, values_fill=0)), row.names=FALSE)
write_csv(Q %>% select(locus, label, csp, vote, gap, S1, S2, S3),
          "DR/out/DR05e_quartets.csv")

hr("3. ADDITIVITY")
cat("  Buneman is exact only for additive distances. Near-ties carry no\n")
cat("  topological information, so the test is run at several gap cuts.\n\n")
print(as.data.frame(Q %>% group_by(label) %>%
  summarise(quartets=n(), median_gap=round(median(gap),4),
            f01=round(mean(gap>0.01),3), f02=round(mean(gap>0.02),3),
            f05=round(mean(gap>0.05),3), .groups="drop")), row.names=FALSE)
p2 <- ggplot(Q, aes(gap, fill=label)) +
  geom_histogram(bins=60, alpha=0.6, position="identity", colour=NA) +
  scale_fill_manual(values=c(A="#1D9E75", B="#D85A30")) +
  xlim(0,0.3) + labs(title="DR05e - quartet resolution",
    subtitle="mass at zero = near-ties, no topological information",
    x="gap = (S_mid - S_min)/mean(S)", y="quartets") + theme_minimal(10)
suppressWarnings(ggsave("DR/fig/DR05e_2_gap.pdf", p2, width=8, height=5))

hr("4. THE TEST — locus is the unit, not the quartet")
cat("  Each locus contributes ONE vote: the majority across its quartets.\n")
cat("  Ties are dropped. Binomial against 1/3, the null when the three\n")
cat("  resolutions are equally likely.\n\n")
for (GAPCUT in c(0, 0.01, 0.02, 0.05)) {
  q <- Q %>% filter(gap > GAPCUT)
  L <- q %>% count(locus, vote) %>% group_by(locus) %>%
    slice_max(n, n=1, with_ties=TRUE) %>%
    summarise(vote = if (n() == 1) vote else NA_character_,
              nq = sum(n), .groups="drop") %>% filter(!is.na(vote))
  tb <- table(factor(L$vote, levels=c("regia+Dionaea","regia+core","neither")))
  k <- tb[["regia+Dionaea"]]; n <- sum(tb)
  bt <- binom.test(k, n, 1/3)
  cat(sprintf("  gap > %.2f : loci %4d | regia+Dionaea %4d (%.3f) | regia+core %4d | neither %4d\n",
              GAPCUT, n, k, k/n, tb[["regia+core"]], tb[["neither"]]))
  cat(sprintf("              binomial vs 1/3: p = %.3g | 95%% CI [%.3f, %.3f]\n",
              bt$p.value, bt$conf.int[1], bt$conf.int[2]))
  if (GAPCUT == 0.02) LOCI <- L
}
cat("\n  Also per subgenome separately, since A and B are independent evidence:\n")
for (SUB in c("A","B")) {
  q <- Q %>% filter(gap > 0.02, label == SUB)
  L <- q %>% count(locus, vote) %>% group_by(locus) %>%
    slice_max(n, n=1, with_ties=TRUE) %>%
    summarise(vote = if (n() == 1) vote else NA_character_, .groups="drop") %>%
    filter(!is.na(vote))
  tb <- table(factor(L$vote, levels=c("regia+Dionaea","regia+core","neither")))
  bt <- binom.test(tb[["regia+Dionaea"]], sum(tb), 1/3)
  cat(sprintf("    %s : %d loci | regia+Dionaea %.3f | p = %.3g\n",
              SUB, sum(tb), tb[["regia+Dionaea"]]/sum(tb), bt$p.value))
}
write_csv(LOCI, "DR/out/DR05e_loci.csv")

hr("5. CONSISTENCY across core species")
cat("  If regia really is sister to Dionaea, EVERY core species should say so.\n")
cat("  If only one or two do, something species-specific is driving it.\n\n")
print(as.data.frame(Q %>% filter(gap > 0.02) %>%
  group_by(csp, label) %>%
  summarise(quartets=n(), regia_Dionaea=round(mean(vote=="regia+Dionaea"),3),
            regia_core=round(mean(vote=="regia+core"),3),
            neither=round(mean(vote=="neither"),3), .groups="drop")),
  row.names=FALSE)
p1 <- ggplot(Q %>% filter(gap > 0.02), aes(csp, fill=vote)) +
  geom_bar(position="fill") +
  geom_hline(yintercept=1/3, colour="grey30", linetype="dashed") +
  facet_wrap(~label) +
  scale_fill_manual(values=c("regia+Dionaea"="#D85A30",
                             "regia+core"="#1D9E75", "neither"="grey65")) +
  labs(title="DR05e - rate-free quartet votes on regia's position",
       subtitle="dashed 1/3 = null | per-quartet view; the TEST is on locus majorities",
       x="core species in the quartet", y="fraction of quartets", fill=NULL) +
  theme_minimal(10) + theme(legend.position="top",
                            axis.text.x=element_text(angle=30, hjust=1))
suppressWarnings(ggsave("DR/fig/DR05e_1_votes.pdf", p1, width=9, height=5))

hr("VERDICT")
cat("  The four-point cannot be fooled by regia being slow: each tip enters\n")
cat("  each of the three sums exactly once, so rate cancels from the ranking.\n\n")
cat("  regia+Dionaea well above 1/3 in BOTH subgenomes -> rate is excluded and\n")
cat("    regia is sister to Dionaea. Drosera would be PARAPHYLETIC.\n")
cat("  regia+core above 1/3 -> the concatenated and ASTRAL results are a rate\n")
cat("    artefact after all, and regia is an early-branching Drosera.\n")
cat("  all three near 1/3 -> no rate-free signal; the node is unresolvable and\n")
cat("    should be reported as such.\n\n")
cat("  NOTE BEFORE REPORTING: Drosera regia is conventionally the earliest-\n")
cat("  diverging Drosera, in its own subgenus. Sister to Dionaea would make\n")
cat("  Drosera paraphyletic -- a much stronger claim. Check the literature.\n")
