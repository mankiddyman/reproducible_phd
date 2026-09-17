#!/usr/bin/env Rscript
# ============================================================================
# DR08 — IS THE A-DUPLICATION IN CORE DROSERA AUTO OR ALLO?
#
# DR07 established the A-duplication is SHARED across binata, paradoxa,
# scorpioides and capensis (0.000-0.099 independent, p to 6e-22), and NOT
# shared with regia (pooled 0.385, CI [0.325,0.448], at the null).
# So those four inherited two A sublineages from a common ancestor.
# This asks HOW that second A arose.
#
# THE QUARTET  (copy1, copy2, Dionaea_A, Nepenthes)
#   (c1,c2)|(DioA,Nep)   the copies are SISTERS, nothing between them
#                        -> AUTOpolyploid: one lineage doubled
#   (c1,DioA)|(c2,Nep)   Dionaea branches BETWEEN the copies
#   (c1,Nep)|(c2,DioA)   -> ALLOpolyploid: the copies descend from lineages
#                           that split BEFORE the Drosera/Dionaea divergence
#
#   NULL: 1/3 auto, 2/3 allo, since two of three resolutions mean allo.
#   Copy 1 and copy 2 are arbitrary labels WITHIN a locus, so no cross-locus
#   phasing is needed. That is why this is cheap and phasing is not.
#
# THIS IS SCRIPT 17 ONE LEVEL DOWN. There it asked whether Dionaea's A and B
# are sisters -- they are not, 17.3% vs a 1/3 null, p = 4.1e-19, which is the
# evidence for the original allopolyploidy. Same machinery, same logic.
#
# WHY IT MATTERS FOR DATING
#   AUTO: the two copies were IDENTICAL at the doubling, so their divergence
#         dates the event directly.
#   ALLO: they were already diverged when they met, so their divergence only
#         gives an UPPER BOUND and the hybridisation must be bracketed.
#
# CONTROLS
#   1. regia, which DR07 says does NOT share the duplication. Its answer
#      should differ, or the test is not discriminating.
#   2. Dionaea's own A vs B through the same code, which must reproduce
#      script 17's ~17% sisters. If it does not, the implementation is wrong.
#
# IN   DR/out/DR02_gene_delta.csv, DR/out/DR02_segments.csv,
#      DR/out/pairwise_ks.csv, DR/locus_meta.tsv, fractionation_by_chrpair.csv
# OUT  DR/out/DR08_alloauto.csv
# FIG  DR/fig/DR08_1_alloauto.pdf, DR08_2_depths.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))
COLLAPSE <- 0.25; GAPCUT <- 0.02
CORE <- c("binata","paradoxa","scorpioides","capensis")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

hr("1. SETUP")
K <- read_csv("DR/out/pairwise_ks.csv", show_col_types=FALSE) %>%
     filter(!is.na(dS), dS>=0, dS<5, codons>=100)
DD <- bind_rows(K %>% transmute(locus=anchor, a=seq1, b=seq2, d=dS),
                K %>% transmute(locus=anchor, a=seq2, b=seq1, d=dS))
dmap <- setNames(DD$d, paste(DD$locus, DD$a, DD$b))
gd <- function(l,x,y) ifelse(x==y, 0, unname(dmap[paste(l,x,y)]))
LOC <- read_tsv("DR/locus_meta.tsv", show_col_types=FALSE) %>%
       group_by(tip) %>% slice(1) %>% ungroup()
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
base <- dio %>% inner_join(nep, by="locus")
cat(sprintf("  loci with a Dionaea homeolog pair AND Nepenthes: %d\n", nrow(base)))

G  <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE)
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE)
lab <- SS %>% filter(label %in% c("A","B")) %>%
       select(genome, region, chr, mb_lo, mb_hi, label)
GL <- G %>% mutate(mb=mid/1e6) %>%
  inner_join(lab, by=c("genome","region","chr"), relationship="many-to-many") %>%
  filter(mb>=mb_lo, mb<=mb_hi) %>% group_by(tip) %>% slice(1) %>% ungroup() %>%
  mutate(sp=sub("Drosera_","",genome))
cat("  collapsing recent duplicates at dS <", COLLAPSE, "\n")
GL <- GL %>% group_by(locus, sp, label) %>% group_modify(function(d, key) {
  if (nrow(d) <= 1) return(d)
  n <- nrow(d); cl <- seq_len(n)
  for (i in 1:(n-1)) for (j in (i+1):n) {
    dd <- gd(key$locus, d$tip[i], d$tip[j])
    if (!is.na(dd) && dd < COLLAPSE) cl[cl==cl[j]] <- cl[i]
  }
  d[!duplicated(cl), ]
}) %>% ungroup()

# the shared quartet engine
qtest <- function(l, c1, c2, out, ogr) {
  S <- c(gd(l,c1,c2) + gd(l,out,ogr),
         gd(l,c1,out) + gd(l,c2,ogr),
         gd(l,c1,ogr) + gd(l,c2,out))
  if (any(is.na(S))) return(c(NA_character_, NA_real_))
  o <- sort(S); g <- (o[2]-o[1])/mean(S)
  c(c("sisters","split","split")[which.min(S)], g)
}
run <- function(tips_fun, label_txt) {
  r <- bind_rows(lapply(base$locus, function(l) {
    tp <- tips_fun(l); if (is.null(tp)) return(NULL)
    v <- qtest(l, tp[1], tp[2], base$DioA[base$locus==l], base$Nep[base$locus==l])
    if (is.na(v[1])) return(NULL)
    tibble(locus=l, res=v[1], gap=as.numeric(v[2]),
           d_copies=gd(l, tp[1], tp[2]))
  }))
  if (!nrow(r)) return(NULL)
  rr <- r %>% filter(gap > GAPCUT)
  k <- sum(rr$res=="sisters"); n <- nrow(rr)
  bt <- binom.test(k, n, 1/3)
  tibble(test=label_txt, n=n, sisters=round(k/n,3),
         split=round(1-k/n,3),
         CI=sprintf("[%.2f,%.2f]", bt$conf.int[1], bt$conf.int[2]),
         p_vs_third=signif(bt$p.value,3),
         med_d_copies=round(median(rr$d_copies, na.rm=TRUE),3),
         verdict=ifelse(bt$p.value > 0.05, "ambiguous",
                 ifelse(k/n > 1/3, "AUTO (copies are sisters)",
                        "ALLO (Dionaea branches between)")))
}

hr("2. CONTROL — reproduce script 17 on Dionaea's own A vs B")
cat("  Script 17 found Dionaea's copies are sisters at 17.3% vs a 1/3 null,\n")
cat("  p = 4.1e-19. If this code does not reproduce that, it is wrong.\n\n")
ctl <- run(function(l) {
  a <- base$DioA[base$locus==l]; b <- base$DioB[base$locus==l]
  if (!length(a) || !length(b)) NULL else c(a, b)
}, "Dionaea A vs B  [must give ~0.17]")
# NOTE: for this control the 'outgroup' slot cannot also be DioA -- rebuild
ctl <- local({
  r <- bind_rows(lapply(base$locus, function(l) {
    a <- base$DioA[base$locus==l]; b <- base$DioB[base$locus==l]
    n_ <- base$Nep[base$locus==l]
    dr <- GL$tip[GL$locus==l & GL$sp=="binata" & GL$label=="A"]
    if (!length(dr)) return(NULL)
    v <- qtest(l, a, b, dr[1], n_)
    if (is.na(v[1])) return(NULL)
    tibble(res=v[1], gap=as.numeric(v[2]))
  }))
  rr <- r %>% filter(gap > GAPCUT)
  k <- sum(rr$res=="sisters"); n <- nrow(rr)
  bt <- binom.test(k, n, 1/3)
  tibble(test="Dionaea A vs B, binata_A as the third tip", n=n,
         sisters=round(k/n,3), split=round(1-k/n,3),
         CI=sprintf("[%.2f,%.2f]", bt$conf.int[1], bt$conf.int[2]),
         p_vs_third=signif(bt$p.value,3), med_d_copies=NA_real_,
         verdict=ifelse(k/n < 1/3, "ALLO -- matches script 17", "MISMATCH"))
})
print(as.data.frame(ctl), row.names=FALSE)

hr("3. THE TEST — core Drosera's two A copies")
cat("  Dionaea_A is the third tip; Nepenthes the outgroup.\n\n")
res <- bind_rows(lapply(c(CORE, "regia"), function(s) {
  run(function(l) {
    t <- GL$tip[GL$locus==l & GL$sp==s & GL$label=="A"]
    if (length(t) != 2) NULL else t
  }, s)
}))
print(as.data.frame(res), row.names=FALSE)
cat("\n  sisters > 1/3 -> AUTO: one lineage doubled, copies identical at the event.\n")
cat("  sisters < 1/3 -> ALLO: Dionaea branches between them, so the two A\n")
cat("    lineages had already split before Drosera and Dionaea diverged.\n")
cat("    That is a SECOND HYBRIDISATION, not a whole-genome doubling.\n")
write_csv(bind_rows(ctl, res), "DR/out/DR08_alloauto.csv")

hr("4. DEPTH CHECK — does it agree with the topology?")
cat("  If ALLO, the two A copies should be nearly as divergent as A is from\n")
cat("  Dionaea's A. If AUTO, much shallower.\n\n")
dep <- bind_rows(lapply(c(CORE, "regia"), function(s) {
  d <- bind_rows(lapply(base$locus, function(l) {
    t <- GL$tip[GL$locus==l & GL$sp==s & GL$label=="A"]
    if (length(t) != 2) return(NULL)
    a <- base$DioA[base$locus==l]
    tibble(sp=s, d_within=gd(l,t[1],t[2]),
           d_to_dio=mean(c(gd(l,t[1],a), gd(l,t[2],a))))
  }))
  d %>% filter(!is.na(d_within), !is.na(d_to_dio), d_to_dio>0) %>%
    summarise(sp=s[1], loci=n(),
              med_within=round(median(d_within),3),
              med_to_Dionaea=round(median(d_to_dio),3),
              ratio=round(median(d_within/d_to_dio),3))
}))
print(as.data.frame(dep), row.names=FALSE)
cat("\n  ratio near 1 -> the copies are as distant from each other as from\n")
cat("  Dionaea, i.e. ALLO. Ratio well below 1 -> AUTO.\n")

p1 <- ggplot(res, aes(reorder(test, sisters), sisters,
                      fill=test=="regia")) +
  geom_hline(yintercept=1/3, linetype="dashed", colour="grey40") +
  geom_col(width=0.6) +
  geom_text(aes(label=sprintf("%.3f  (n=%d)", sisters, n)), hjust=-0.12, size=3) +
  annotate("text", x=0.6, y=0.345, hjust=0, size=3, colour="grey35",
           label="1/3 = null") +
  scale_fill_manual(values=c("FALSE"="#1D9E75","TRUE"="grey70"), guide="none") +
  coord_flip(ylim=c(0,0.8)) +
  labs(title="DR08 - is the A-duplication a WGD or a second hybridisation?",
       subtitle=paste0("Quartet (copy1, copy2, Dionaea_A, Nepenthes). ",
         "Above 1/3 = copies are sisters = AUTOpolyploid WGD.\n",
         "Below = Dionaea branches between them = ALLOpolyploid, ",
         "a second hybridisation."),
       x=NULL, y="fraction of loci where the two A copies are sisters") +
  theme_minimal(11) + theme(plot.subtitle=element_text(size=8.5, colour="grey30"))
suppressWarnings(ggsave("DR/fig/DR08_1_alloauto.pdf", p1, width=9, height=5))

hr("VERDICT")
cat("  AUTO -> a straightforward doubling. The copies' divergence DATES the\n")
cat("    event, which is the easy case for DR09.\n")
cat("  ALLO -> a SECOND HYBRIDISATION in core Drosera, with an A-lineage that\n")
cat("    had already diverged. Then the copy divergence is only an UPPER\n")
cat("    BOUND and the event must be bracketed, exactly as for the original\n")
cat("    A/B hybridisation. It also makes phasing worth doing, because you\n")
cat("    would want to know whether the added genome sits in coherent blocks.\n")
cat("\n  regia is the discriminating control: DR07 says it does NOT share this\n")
cat("  duplication, so its answer should differ from the core four.\n")
