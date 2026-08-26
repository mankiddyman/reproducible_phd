#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); BAND <- 0.08
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
CORE <- c("Drosera_binata","Drosera_paradoxa","Drosera_scorpioides")
REFS <- c(CORE, "Drosera_regia")            # capensis dropped as a reference

k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS>0, dS<5, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
key <- tm %>% select(anchor=nep_gene, k1=tip, genome, chr) %>% distinct() %>%
       mutate(region=sub("-.*$","",anchor))
DD <- bind_rows(k %>% transmute(anchor, ga=seq1, gb=seq2, d=dS),
                k %>% transmute(anchor, ga=seq2, gb=seq1, d=dS))
seg <- read_csv("trees/qc/segment_delta_calibrated.csv", show_col_types=FALSE) %>%
  mutate(side=case_when(med_delta < -BAND ~ "A", med_delta > BAND ~ "B", TRUE ~ "amb")) %>%
  select(genome, region, chr, side)

assign_one <- function(FO, RE) {
  f <- key %>% filter(genome==FO) %>% select(anchor, region, fg=k1, fchr=chr)
  r <- key %>% filter(genome==RE) %>% group_by(anchor) %>%
       filter(n_distinct(chr) >= 2) %>% ungroup() %>% select(anchor, rg=k1, rchr=chr)
  f %>% inner_join(r, by="anchor", relationship="many-to-many") %>%
    inner_join(DD %>% select(anchor, fg=ga, rg=gb, d), by=c("anchor","fg","rg")) %>%
    group_by(anchor, region, fg, fchr) %>% arrange(d, .by_group=TRUE) %>%
    summarise(best_rchr=rchr[1], margin=if (n()>1) (d[2]-d[1])/d[1] else NA_real_,
              .groups="drop") %>% mutate(focal=FO, ref=RE)
}
build <- function(refset) {
  pr <- expand.grid(f=CORE, r=refset, stringsAsFactors=FALSE) %>% filter(f != r)
  bind_rows(lapply(seq_len(nrow(pr)), function(i) assign_one(pr$f[i], pr$r[i]))) %>%
    filter(!is.na(margin), margin > 0.5) %>%
    left_join(seg %>% rename(ref=genome, best_rchr=chr, rside=side),
              by=c("ref","region","best_rchr")) %>%
    filter(!is.na(rside), rside != "amb")
}
consensus <- function(AS) AS %>% group_by(focal, anchor, fg) %>%
  summarise(side=names(sort(table(rside), decreasing=TRUE))[1],
            unanimous=n_distinct(rside)==1, nref=n_distinct(ref), .groups="drop")

hr("1. reference-set comparison")
for (nm in c("core+regia (capensis dropped)","core only","all five (previous run)")) {
  rs <- switch(nm, "core only"=CORE, "all five (previous run)"=c(CORE,"Drosera_regia","Drosera_capensis"), REFS)
  cn <- consensus(build(rs))
  cat(sprintf("\n  %s\n", nm))
  print(as.data.frame(cn %>% group_by(focal) %>%
    summarise(copies=n(), A=sum(side=="A"), B=sum(side=="B"),
              ratio=round(sum(side=="A")/pmax(1,sum(side=="B")),2), .groups="drop")),
    row.names=FALSE)
}

AS <- build(REFS); cons <- consensus(AS)

hr("2. leave-one-reference-out stability")
print(as.data.frame(bind_rows(lapply(REFS, function(drop) {
  cn <- consensus(build(setdiff(REFS, drop)))
  cn %>% group_by(focal) %>%
    summarise(ratio=round(sum(side=="A")/pmax(1,sum(side=="B")),2), .groups="drop") %>%
    mutate(dropped=sub("Drosera_","",drop))
})) %>% pivot_wider(names_from=dropped, values_from=ratio)), row.names=FALSE)
cat("\n  Stable across drops = no single reference is driving the ratio.\n")

hr("3. THE A=0 TEST — is there exactly ONE B subgenome?")
cat("  AAB has one B, so a 2-copy locus can never be B+B: A=0 impossible.\n")
cat("  ABB has two B, so A=2 becomes impossible instead.\n")
cat("  Free binomial with the observed p_A predicts both at (1-p)^2 and p^2.\n\n")
loc <- cons %>% group_by(focal, anchor) %>%
       summarise(n=n(), A=sum(side=="A"), .groups="drop")
two <- loc %>% filter(n==2)
print(as.data.frame(bind_rows(lapply(CORE, function(f) {
  d <- two %>% filter(focal==f); cn <- cons %>% filter(focal==f)
  p <- mean(cn$side=="A"); n <- nrow(d)
  tibble(focal=sub("Drosera_","",f), loci=n, p_A=round(p,3),
         A0=sum(d$A==0), A0_exp=round((1-p)^2*n,1),
         p_A0=signif(binom.test(sum(d$A==0), n, (1-p)^2)$p.value,3),
         A2=sum(d$A==2), A2_exp=round(p^2*n,1),
         p_A2=signif(binom.test(sum(d$A==2), n, p^2)$p.value,3))
}))), row.names=FALSE)
cat("\n  A=0 far below expectation with A=2 at or above it => exactly one B.\n")

hr("4. 2-copy composition against the strict-AAB prediction")
cat("  Strict AAB with random retention: A=1 twice as common as A=2.\n\n")
print(as.data.frame(bind_rows(lapply(CORE, function(f) {
  d <- two %>% filter(focal==f); a1 <- sum(d$A==1); a2 <- sum(d$A==2)
  tibble(focal=sub("Drosera_","",f), A0=sum(d$A==0), A1=a1, A2=a2,
         obs_ratio=round(a1/pmax(1,a2),2),
         p_vs_2to1=signif(binom.test(a1, a1+a2, 2/3)$p.value,3))
}))), row.names=FALSE)
cat("\n  Below 2.0 = more all-A loci than strict AAB allows. Consistent with\n")
cat("  differential fractionation, or with more than two A subgenomes.\n")

hr("5. 3-copy composition")
print(as.data.frame(loc %>% filter(n==3) %>% group_by(focal, A) %>%
  summarise(loci=n(), .groups="drop") %>%
  pivot_wider(names_from=A, values_from=loci, values_fill=0, names_prefix="A=")),
  row.names=FALSE)
cat("\n  Strict AAB predicts every 3-copy locus at A=2.\n")

hr("6. how accurate is a single side call?")
ag <- AS %>% group_by(focal, anchor, fg) %>% filter(n_distinct(ref)>=2) %>%
  summarise(agree=n_distinct(rside)==1, .groups="drop")
a <- mean(ag$agree); e <- (2 - sqrt(4 - 8*(1-a)))/4
cat(sprintf("  two-reference agreement: %.3f  => implied per-call error ~%.3f\n", a, e))
cat(sprintf("  under strict AAB a 3-copy locus reads 2A:1B with prob (1-e)^3 = %.2f\n",
            (1-e)^3))
cat(sprintf("  observed: %s\n", paste(sprintf("%s %.2f", sub("Drosera_","",CORE),
  vapply(CORE, function(f){d <- loc %>% filter(focal==f, n==3); mean(d$A==2)}, numeric(1))),
  collapse="  ")))
cat("\n  CAVEAT: both references inherit their A/B labels from the same Dionaea\n")
cat("  delta axis, so their errors are correlated and this understates e.\n")
write_csv(cons, "trees/qc/consensus_norefcap.csv")
