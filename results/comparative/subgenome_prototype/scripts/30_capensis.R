#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); BAND <- 0.08
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
CAP <- "Drosera_capensis"
REFS <- c("Drosera_binata","Drosera_paradoxa","Drosera_scorpioides","Drosera_regia")

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

hr("0. capensis own-copy dS — where to cut the recent duplication")
wi <- k %>% filter(sp1==CAP, sp2==CAP) %>%
  left_join(key %>% select(anchor, seq1=k1, c1=chr), by=c("anchor","seq1")) %>%
  left_join(key %>% select(anchor, seq2=k1, c2=chr), by=c("anchor","seq2")) %>%
  filter(!is.na(c1), !is.na(c2), c1 != c2)
cat(sprintf("  different-chromosome pairs: %d\n", nrow(wi)))
print(round(quantile(wi$dS, c(.05,.1,.2,.25,.3,.4,.5,.75,.9)), 3))
cat("\n  FIG46 gave capensis mode_min 0.092 / med_max 0.723; 23 sec3 gave\n")
cat("  same-B 0.106, same-A 0.640, cross-A/B 1.416. Cut between 0.1 and 0.7.\n")

comps <- function(ids, e, TT) {
  lab <- setNames(seq_along(ids), ids); ee <- e[e$d < TT, , drop=FALSE]
  if (nrow(ee)) for (it in seq_along(ids)) for (j in seq_len(nrow(ee))) {
    a <- lab[[ee$g1[j]]]; b <- lab[[ee$g2[j]]]
    if (a != b) { m <- min(a,b); lab[lab %in% c(a,b)] <- m } }
  lab
}
cluster <- function(SP, TT) {
  PW <- k %>% filter(sp1==SP, sp2==SP) %>% transmute(anchor, g1=seq1, g2=seq2, d=dS)
  key %>% filter(genome==SP) %>% group_by(anchor) %>% filter(n()>=2) %>%
    group_modify(function(d, kk){
      ids <- d$k1
      e <- PW[PW$anchor==kk$anchor & PW$g1 %in% ids & PW$g2 %in% ids, , drop=FALSE]
      tibble(k1=ids, chr=d$chr, region=d$region, cl=unname(comps(ids, e, TT)))
    }) %>% ungroup()
}

hr("1. collapse sweep — clusters per locus")
cat("  capensis is dodecaploid; v2 8.3 gives median 5 major tracts (3-6).\n")
cat("  Too low a threshold leaves recent duplicates as separate clusters;\n")
cat("  too high merges ancestral subgenomes.\n\n")
for (TT in c(0.15, 0.25, 0.35, 0.5, 0.7)) {
  cl <- cluster(CAP, TT)
  nl <- cl %>% group_by(anchor) %>% summarise(ncl=n_distinct(cl), ncop=n(), .groups="drop")
  cat(sprintf("  T=%.2f | loci %d | mean copies %.2f | mean clusters %.2f | ",
              TT, nrow(nl), mean(nl$ncop), mean(nl$ncl)))
  print(table(nl$ncl))
}

## assign a CLUSTER to A or B by its nearest reference copy (mean over cluster members)
assign_cl <- function(cl, RE) {
  r <- key %>% filter(genome==RE) %>% group_by(anchor) %>%
       filter(n_distinct(chr) >= 2) %>% ungroup() %>% select(anchor, rg=k1, rchr=chr)
  cl %>% inner_join(r, by="anchor", relationship="many-to-many") %>%
    inner_join(DD %>% select(anchor, k1=ga, rg=gb, d), by=c("anchor","k1","rg")) %>%
    group_by(anchor, region, cl, rchr) %>%
    summarise(dmean=mean(d), .groups="drop") %>%
    group_by(anchor, region, cl) %>% arrange(dmean, .by_group=TRUE) %>%
    summarise(best_rchr=rchr[1],
              margin=if (n()>1) (dmean[2]-dmean[1])/dmean[1] else NA_real_,
              .groups="drop") %>% mutate(ref=RE)
}
run <- function(SP, TT) {
  cl <- cluster(SP, TT)
  bind_rows(lapply(REFS[REFS != SP], function(R) assign_cl(cl, R))) %>%
    filter(!is.na(margin), margin > 0.5) %>%
    left_join(seg %>% rename(ref=genome, best_rchr=chr, rside=side),
              by=c("ref","region","best_rchr")) %>%
    filter(!is.na(rside), rside != "amb") %>%
    group_by(anchor, region, cl) %>%
    summarise(side=names(sort(table(rside), decreasing=TRUE))[1],
              nref=n_distinct(ref), unanimous=n_distinct(rside)==1, .groups="drop") %>%
    mutate(focal=SP, T=TT)
}

hr("2. VALIDATION — same procedure on species with a known 2:1 answer")
cat("  Core species have no recent duplication, so collapsing should not\n")
cat("  change their answer. If it does, the threshold is over-merging.\n\n")
print(as.data.frame(bind_rows(lapply(c("Drosera_binata","Drosera_paradoxa",
                                       "Drosera_scorpioides"), function(s)
  bind_rows(lapply(c(0.25, 0.35), function(TT) run(s, TT))))) %>%
  group_by(focal, T) %>%
  summarise(clusters=n(), A=sum(side=="A"), B=sum(side=="B"),
            ratio=round(sum(side=="A")/pmax(1,sum(side=="B")),2), .groups="drop")),
  row.names=FALSE)
cat("\n  Expected: binata 2.04, paradoxa 2.03, scorpioides 2.22 (script 29).\n")

hr("3. CAPENSIS CONSTITUTION")
CP <- bind_rows(lapply(c(0.15, 0.25, 0.35, 0.5), function(TT) run(CAP, TT)))
print(as.data.frame(CP %>% group_by(T) %>%
  summarise(loci=n_distinct(anchor), clusters=n(),
            A=sum(side=="A"), B=sum(side=="B"),
            ratio=round(sum(side=="A")/pmax(1,sum(side=="B")),2),
            unanim=round(mean(unanimous),3), .groups="drop")), row.names=FALSE)
cat("\n  Doubled AAB (AAAABB) predicts 2.0. Stable across T = real.\n")

hr("4. per-locus cluster composition, T = 0.35")
loc <- CP %>% filter(T==0.35) %>% group_by(anchor) %>%
  summarise(n=n(), A=sum(side=="A"), B=sum(side=="B"), .groups="drop")
print(as.data.frame(loc %>% group_by(n, A) %>% summarise(loci=n(), .groups="drop") %>%
  pivot_wider(names_from=A, values_from=loci, values_fill=0, names_prefix="A=")),
  row.names=FALSE)
p <- sum(loc$A)/sum(loc$n)
cat(sprintf("\n  p_A = %.3f\n", p))
for (nn in 2:4) {
  d <- loc %>% filter(n==nn); if (nrow(d) < 15) next
  cat(sprintf("  %d-cluster loci: %d | A=0 obs %d exp %.1f | p = %.3g\n", nn, nrow(d),
    sum(d$A==0), (1-p)^nn*nrow(d),
    binom.test(sum(d$A==0), nrow(d), (1-p)^nn)$p.value))
}
cat("\n  A=0 depletion measures how many B subgenomes there are, as in 29 sec3.\n")

hr("5. was the recent duplication A-biased or B-biased?")
cat("  Cluster size = how many copies collapsed into one ancestral lineage.\n\n")
cl35 <- cluster(CAP, 0.35) %>% count(anchor, cl, name="size")
print(as.data.frame(CP %>% filter(T==0.35) %>%
  inner_join(cl35, by=c("anchor","cl")) %>% group_by(side) %>%
  summarise(clusters=n(), mean_size=round(mean(size),3),
            frac_dup=round(mean(size>1),3), .groups="drop")), row.names=FALSE)
cat("\n  Similar frac_dup on both sides = the extra WGD doubled the whole genome.\n")
cat("  Strongly unequal = it hit one ancestry preferentially.\n")
write_csv(CP, "trees/qc/capensis_constitution.csv")
