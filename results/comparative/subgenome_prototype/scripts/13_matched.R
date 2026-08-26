#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); NPERM <- 999
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"

k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS > 0, dS < 5, codons >= 100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
cn <- tm %>% count(nep_gene, genome, name="copies")
sc <- cn %>% filter(copies == 1)

SP <- k %>% filter(sp1!=sp2, sp1!=NEP, sp2!=NEP) %>%
  semi_join(sc, by=c("anchor"="nep_gene","sp1"="genome")) %>%
  semi_join(sc, by=c("anchor"="nep_gene","sp2"="genome")) %>%
  mutate(a=pmin(sp1,sp2), b=pmax(sp1,sp2)) %>% select(a,b,dS)
WI <- k %>% filter(sp1==sp2, sp1!=NEP) %>% select(sp=sp1, dS)

hr("1. matched-pair test — loci where BOTH species have exactly 2 copies")
cat("  Optimal 2x2 matching. d2 = larger of the two matched distances.\n")
cat("  Two lineages shared -> d2 ~ speciation.  One lineage duplicated -> d2 ~ ancient.\n\n")
X <- k %>% filter(sp1!=sp2, sp1!=NEP, sp2!=NEP) %>%
  mutate(a=pmin(sp1,sp2), b=pmax(sp1,sp2),
         ga=ifelse(sp1==a, seq1, seq2), gb=ifelse(sp1==a, seq2, seq1)) %>%
  inner_join(cn %>% rename(a=genome, ka=copies), by=c("anchor"="nep_gene","a")) %>%
  inner_join(cn %>% rename(b=genome, kb=copies), by=c("anchor"="nep_gene","b")) %>%
  filter(ka==2, kb==2)

md <- X %>% group_by(a, b, anchor) %>% filter(n()==4) %>%
  group_modify(function(d, key) {
    ia <- sort(unique(d$ga)); ib <- sort(unique(d$gb))
    M <- matrix(NA_real_, 2, 2)
    for (r in 1:2) for (cc in 1:2) {
      v <- d$dS[d$ga==ia[r] & d$gb==ib[cc]]; if (length(v)) M[r,cc] <- v[1] }
    if (any(is.na(M))) return(tibble())
    s1 <- M[1,1]+M[2,2]; s2 <- M[1,2]+M[2,1]
    p <- if (s1 <= s2) c(M[1,1],M[2,2]) else c(M[1,2],M[2,1])
    tibble(d1=min(p), d2=max(p), unmatched=mean(setdiff(as.vector(M), p)))
  }) %>% ungroup()

set.seed(1)
out <- bind_rows(lapply(split(md, paste(md$a, md$b)), function(d) {
  pool <- SP$dS[SP$a==d$a[1] & SP$b==d$b[1]]
  if (length(pool) < 30 || nrow(d) < 20) return(NULL)
  spec <- median(pool)
  nul <- vapply(seq_len(NPERM), function(i)
    median(vapply(seq_len(nrow(d)), function(j) max(sample(pool,2,replace=TRUE)), numeric(1))),
    numeric(1))
  wi <- median(c(WI$dS[WI$sp==d$a[1]], WI$dS[WI$sp==d$b[1]]))
  tibble(a=d$a[1], b=d$b[1], loci=nrow(d), spec=round(spec,3),
         d1=round(median(d$d1),3), d2=round(median(d$d2),3),
         d2_null=round(median(nul),3),
         null_hi=round(quantile(nul,.975),3),
         within_sp=round(wi,3),
         d2_vs_null=round(median(d$d2)/median(nul),3),
         verdict=ifelse(median(d$d2) <= quantile(nul,.975), "two lineages", "check"))
})) %>% arrange(d2_vs_null)
print(as.data.frame(out), row.names=FALSE)
cat("\n  d2 near d2_null  = both copies have an ortholog partner = TWO lineages.\n")
cat("  d2 near within_sp = only one lineage, duplicated.\n")

hr("2. binata specifically")
bn <- out %>% filter(a=="Drosera_binata" | b=="Drosera_binata")
print(as.data.frame(bn %>% select(a,b,loci,spec,d2,d2_null,within_sp,d2_vs_null,verdict)),
      row.names=FALSE)
cat("\n  binata within-species dS mode is 0.820; speciation to paradoxa 0.499.\n")
cat("  If d2 tracks spec rather than 0.82, binata's two copies are A and B.\n")
write_csv(out, "trees/qc/matched_pair_test.csv")

hr("3. capensis confound check — restrict to ancient own-copies")
cat("  Correlated recent duplicates break the independent-draw null in test 12.\n")
cat("  Here: keep only loci where the species' own copies are older than speciation.\n\n")
own <- k %>% filter(sp1==sp2, sp1!=NEP) %>% group_by(anchor, sp=sp1) %>%
       summarise(own_max=max(dS), .groups="drop")
md2 <- md %>%
  inner_join(own %>% rename(a=sp, own_a=own_max), by=c("anchor","a")) %>%
  inner_join(own %>% rename(b=sp, own_b=own_max), by=c("anchor","b"))
print(as.data.frame(bind_rows(lapply(split(md2, paste(md2$a, md2$b)), function(d) {
  pool <- SP$dS[SP$a==d$a[1] & SP$b==d$b[1]]; if (length(pool)<30) return(NULL)
  spec <- median(pool); dd <- d %>% filter(own_a > spec, own_b > spec)
  if (nrow(dd) < 15) return(NULL)
  tibble(a=d$a[1], b=d$b[1], loci_all=nrow(d), loci_anc=nrow(dd),
         d2_all=round(median(d$d2),3), d2_anc=round(median(dd$d2),3),
         spec=round(spec,3), ratio_anc=round(median(dd$d2)/spec,2))
}))), row.names=FALSE)
