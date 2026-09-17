#!/usr/bin/env Rscript
# ============================================================================
# DR11 — WHAT IS THE B-SUBGENOME ASYMMETRY, AND CAN A1/A2 BE PHASED?
#
# PART 1 — THE ASYMMETRY
#   Under PURE incomplete lineage sorting the two minority quartet resolutions
#   must be EQUAL. That is a hard prediction of the coalescent, not an
#   approximation: when lineages fail to coalesce, the order in which they
#   eventually join is random.
#     A subgenome: q2 0.274 vs q3 0.226, z = 1.54  -> symmetric, ILS-consistent
#     B subgenome: q2 0.177 vs q3 0.317, z = 3.68  -> NOT pure ILS
#   ASTRAL does not say which topology q2 and q3 correspond to, so the elevated
#   one is identified here by counting topologies directly in the gene trees.
#
#   The quartet is (regia_B, Dionaea_B, core_B, Nepenthes):
#     T1  (regia,Dionaea)|(core,Nep)   the species tree
#     T2  (regia,core)|(Dionaea,Nep)   regia's B allied with Drosera
#     T3  (regia,Nep)|(Dionaea,core)   Dionaea's B allied with Drosera
#   T2 elevated -> gene flow involving regia's B
#   T3 elevated -> gene flow involving Dionaea's B
#
#   NOTE ON THE D-STATISTIC: T2 vs T3 counted this way IS Patterson's D in
#   topology form -- the same asymmetry ABBA-BABA measures. It is reported as
#   a proportion with a binomial test rather than the block-jackknife Z of the
#   site-pattern version, so treat it as indicative.
#
#   A is run as the control. It should come out symmetric.
#
# PART 2 — CAN A1/A2 BE PHASED?
#   DR07 showed the A-duplication is SHARED among core Drosera and capensis
#   (0.000-0.099 "independent", p to 6e-22) but NOT with regia (0.385 pooled,
#   at the 1/3 null). So the per-locus pairings exist for the core four.
#   What has never been tested is whether those pairings are CONSISTENT along
#   a chromosome -- which is what phasing actually requires.
#
#   The test: assign copies per locus by minimum-distance matching between two
#   species, then ask whether the same physical chromosome pairs with the same
#   chromosome at neighbouring loci. Real subgenomes are inherited as tracts,
#   so consistency should be high WITHIN a chromosome and the assignment should
#   partition by segment. Random interspersion means the pairings are noise.
#
# IN   DR/tree/pro/genetrees_pro_clean.tre, DR/tree/pro/gene2species.map,
#      DR/out/DR02_gene_delta.csv, DR/out/DR02_segments.csv,
#      DR/out/pairwise_ks.csv
# OUT  DR/out/DR11_asymmetry.csv, DR11_phasing.csv
# FIG  DR/fig/DR11_1_asymmetry.pdf, DR11_2_phasing.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(ape); library(dplyr); library(readr); library(tidyr); library(ggplot2)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))
CORE <- c("binata","paradoxa","scorpioides","capensis")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

hr("1. COUNT THE THREE TOPOLOGIES DIRECTLY")
map <- read.table("DR/tree/pro/gene2species.map", sep="\t",
                  col.names=c("gene","unit"), stringsAsFactors=FALSE)
u <- setNames(map$unit, map$gene)
trees <- readLines("DR/tree/pro/genetrees_pro_clean.tre", warn=FALSE)
cat(sprintf("  gene trees: %d\n", length(trees)))

quart <- function(t, a, b, c, d) {
  tips <- c(a,b,c,d)
  if (!all(tips %in% t$tip.label)) return(NA_character_)
  q <- drop.tip(t, setdiff(t$tip.label, tips))
  if (length(q$tip.label) != 4) return(NA_character_)
  sp <- prop.part(unroot(q)); s <- q$tip.label[sp[[length(sp)]]]
  if (length(s) != 2) return(NA_character_)
  paste(sort(s), collapse="|")
}

count_sub <- function(SUB) {
  rg <- paste0("regia_", SUB); di <- paste0("Dionaea_", SUB)
  bind_rows(lapply(CORE, function(cs) {
    co <- paste0(cs, "_", SUB)
    res <- vapply(trees, function(nw) {
      t <- tryCatch(read.tree(text=nw), error=function(e) NULL)
      if (is.null(t)) return(NA_character_)
      # one representative copy per unit, chosen at random per locus
      lab <- t$tip.label; un <- unname(u[lab])
      pick <- function(x) { i <- which(un == x); if (!length(i)) NA else lab[i[1]] }
      a <- pick(rg); b <- pick(di); c <- pick(co); d <- pick("Nepenthes")
      if (any(is.na(c(a,b,c,d)))) return(NA_character_)
      v <- quart(t, a, b, c, d)
      if (is.na(v)) return(NA_character_)
      pr <- strsplit(v, "|", fixed=TRUE)[[1]]
      # a=regia, b=Dionaea, c=core, d=Nepenthes -- these are GENE names
      if (setequal(pr, c(a,b)) || setequal(pr, c(c,d))) "T1_species_tree"
      else if (setequal(pr, c(a,c)) || setequal(pr, c(b,d))) "T2_regia_with_core"
      else "T3_Dionaea_with_core"
    }, character(1))
    res <- res[!is.na(res)]
    tb <- table(factor(res, levels=c("T1_species_tree","T2_regia_with_core",
                                     "T3_Dionaea_with_core")))
    n2 <- tb[["T2_regia_with_core"]]; n3 <- tb[["T3_Dionaea_with_core"]]
    bt <- if (n2+n3 >= 10) binom.test(n2, n2+n3, 0.5) else NULL
    tibble(subgenome=SUB, core=cs, n=sum(tb),
           T1=round(tb[[1]]/sum(tb),3), T2=round(n2/sum(tb),3),
           T3=round(n3/sum(tb),3), n_T2=n2, n_T3=n3,
           D = (n2-n3)/(n2+n3),
           p = if (is.null(bt)) NA_real_ else bt$p.value)
  }))
}
cat("  counting B...\n"); flush.console(); B <- count_sub("B")
cat("  counting A (control)...\n"); flush.console(); A <- count_sub("A")
RES <- bind_rows(B, A) %>% mutate(D=round(D,3), p=signif(p,3))
if (all(RES$T1 == 0) || all(RES$T3 == 1)) {
  cat("\n  *** DEGENERATE: every quartet landed in one category. The\n")
  cat("  classification is broken -- do not read anything below. ***\n")
  quit(status=1)
}
cat("  sanity: T1 should be the MAJORITY (~0.50), matching ASTRAL q1.\n\n")
print(as.data.frame(RES), row.names=FALSE)

hr("2. WHICH ALTERNATIVE IS ELEVATED?")
cat("  D = (T2 - T3)/(T2 + T3). This is Patterson's D in topology form.\n")
cat("    D > 0 -> regia's B allied with core Drosera more often than expected\n")
cat("    D < 0 -> Dionaea's B allied with core Drosera more often\n")
cat("    D ~ 0 -> symmetric, pure ILS\n\n")
pool <- RES %>% group_by(subgenome) %>%
  summarise(quartets=sum(n), T2=sum(n_T2), T3=sum(n_T3), .groups="drop")
pool$D <- round((pool$T2-pool$T3)/(pool$T2+pool$T3), 3)
pool$p <- vapply(seq_len(nrow(pool)), function(i)
  signif(binom.test(pool$T2[i], pool$T2[i]+pool$T3[i], 0.5)$p.value, 3), numeric(1))
print(as.data.frame(pool), row.names=FALSE)
cat("\n  A is the control: it should be near D = 0.\n")
write_csv(RES, "DR/out/DR11_asymmetry.csv")

p1 <- ggplot(RES %>% pivot_longer(c(T1,T2,T3), names_to="topo", values_to="f"),
             aes(core, f, fill=topo)) +
  geom_col(position="stack", width=0.7) +
  geom_hline(yintercept=1/3, linetype="dotted", colour="grey30") +
  facet_wrap(~subgenome, labeller=labeller(subgenome=c(A="A subgenome (control)",
                                                       B="B subgenome"))) +
  scale_fill_manual(values=c(T1_species_tree="#1D9E75",
                             T2_regia_with_core="#D85A30",
                             T3_Dionaea_with_core="grey65"),
                    labels=c("species tree","regia + core","Dionaea + core")) +
  labs(title="DR11 - which disagreement is elevated?",
       subtitle="Pure ILS requires the two minority bars to be EQUAL",
       x="core species in the quartet", y="fraction of gene trees", fill=NULL) +
  theme_minimal(10) + theme(legend.position="top")
suppressWarnings(ggsave("DR/fig/DR11_1_asymmetry.pdf", p1, width=9, height=5))

hr("3. CAN A1/A2 BE PHASED? — consistency along the chromosome")
cat("  DR07 gave per-locus pairings for the core four. Phasing needs those\n")
cat("  pairings to be CONSISTENT along a chromosome, which was never tested.\n\n")
K <- read_csv("DR/out/pairwise_ks.csv", show_col_types=FALSE) %>%
     filter(!is.na(dS), dS>=0, dS<5, codons>=100)
DD <- bind_rows(K %>% transmute(locus=anchor, a=seq1, b=seq2, d=dS),
                K %>% transmute(locus=anchor, a=seq2, b=seq1, d=dS))
dm <- setNames(DD$d, paste(DD$locus, DD$a, DD$b))
gd <- function(l,x,y) ifelse(x==y, 0, unname(dm[paste(l,x,y)]))
G  <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE)
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE)
lab <- SS %>% filter(label=="A") %>% select(genome, region, chr, mb_lo, mb_hi)
GL <- G %>% mutate(mb=mid/1e6) %>%
  inner_join(lab, by=c("genome","region","chr"), relationship="many-to-many") %>%
  filter(mb>=mb_lo, mb<=mb_hi) %>% group_by(tip) %>% slice(1) %>% ungroup() %>%
  mutate(sp=sub("Drosera_","",genome))
two <- GL %>% count(sp, locus, name="k") %>% filter(k==2)
cat("  loci with exactly 2 A copies:\n")
print(as.data.frame(two %>% count(sp, name="loci")), row.names=FALSE)

PH <- bind_rows(lapply(c("binata","paradoxa","scorpioides","regia"), function(s1) {
  s2 <- if (s1=="regia") "binata" else "capensis"
  L <- intersect(two$locus[two$sp==s1], two$locus[two$sp==s2])
  if (length(L) < 20) return(NULL)
  r <- bind_rows(lapply(L, function(l) {
    x <- GL %>% filter(locus==l, sp==s1) %>% arrange(chr, mid)
    y <- GL %>% filter(locus==l, sp==s2) %>% arrange(chr, mid)
    if (nrow(x)!=2 || nrow(y)!=2) return(NULL)
    # minimum-total-distance matching between the two copies
    m11 <- gd(l,x$tip[1],y$tip[1]); m22 <- gd(l,x$tip[2],y$tip[2])
    m12 <- gd(l,x$tip[1],y$tip[2]); m21 <- gd(l,x$tip[2],y$tip[1])
    if (any(is.na(c(m11,m22,m12,m21)))) return(NULL)
    straight <- m11+m22 < m12+m21
    tibble(locus=l,
           chrpair = paste(sort(c(x$chr[1], x$chr[2])), collapse="+"),
           ypair   = paste(sort(c(y$chr[1], y$chr[2])), collapse="+"),
           orient  = if (straight) "straight" else "crossed",
           mid=x$mid[1])
  }))
  if (!nrow(r)) return(NULL)
  cons <- r %>% count(chrpair, ypair, orient) %>%
    group_by(chrpair, ypair) %>%
    summarise(loci=sum(n), top=max(n), consistency=max(n)/sum(n), .groups="drop") %>%
    filter(loci >= 8)
  if (!nrow(cons)) return(NULL)
  tibble(pair=paste(s1,s2,sep="-"), chr_tracks=nrow(cons),
         loci=sum(cons$loci),
         median_consistency=round(median(cons$consistency),3),
         tracks_above_0.8=sum(cons$consistency > 0.8))
}))
if (nrow(PH)) {
  print(as.data.frame(PH), row.names=FALSE)
  cat("\n  consistency = fraction of loci on a chromosome pair whose matching\n")
  cat("  picks the SAME partner chromosome. 0.5 = random. Above ~0.8 means the\n")
  cat("  pairing tracks a real chromosomal unit and phasing is possible.\n")
  write_csv(PH, "DR/out/DR11_phasing.csv")
} else cat("  too few loci with 2 A copies in both species for this test\n")

hr("VERDICT")
cat("  ASYMMETRY: read section 2. If B has D far from 0 while A does not,\n")
cat("  the sign says which lineage's B subgenome is involved.\n\n")
cat("  PHASING: read section 3. median_consistency near 0.5 means the per-locus\n")
cat("  pairings do not track chromosomes and A1/A2 CANNOT be phased. Near 0.9\n")
cat("  means they can, at least for the core four -- regia would still fail\n")
cat("  because DR07 found no shared duplication to align it to.\n")
