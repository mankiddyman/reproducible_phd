#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); BAND <- 0.08
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
CORE <- c("Drosera_binata","Drosera_paradoxa","Drosera_scorpioides")
ALLR <- c(CORE, "Drosera_regia", "Drosera_capensis")

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
cat(sprintf("segment labels: %d  (A %d / B %d / amb %d)\n", nrow(seg),
            sum(seg$side=="A"), sum(seg$side=="B"), sum(seg$side=="amb")))

## nearest reference copy — now returns the GENE, not just its chromosome
assign_one <- function(FO, RE) {
  f <- key %>% filter(genome==FO) %>% select(anchor, region, fg=k1, fchr=chr)
  r <- key %>% filter(genome==RE) %>% group_by(anchor) %>%
       filter(n_distinct(chr) >= 2) %>% ungroup() %>% select(anchor, rg=k1, rchr=chr)
  f %>% inner_join(r, by="anchor", relationship="many-to-many") %>%
    inner_join(DD %>% select(anchor, fg=ga, rg=gb, d), by=c("anchor","fg","rg")) %>%
    group_by(anchor, region, fg, fchr) %>% arrange(d, .by_group=TRUE) %>%
    summarise(best_rg=rg[1], best_rchr=rchr[1], d1=d[1],
              margin=if (n()>1) (d[2]-d[1])/d[1] else NA_real_, .groups="drop") %>%
    mutate(focal=FO, ref=RE)
}
pr <- expand.grid(f=CORE, r=ALLR, stringsAsFactors=FALSE) %>% filter(f != r)
AS <- bind_rows(lapply(seq_len(nrow(pr)), function(i) assign_one(pr$f[i], pr$r[i]))) %>%
      filter(!is.na(margin), margin > 0.5) %>%
      left_join(seg %>% rename(ref=genome, best_rchr=chr, rside=side),
                by=c("ref","region","best_rchr")) %>%
      filter(!is.na(rside), rside != "amb")
cat(sprintf("confident, labelled assignments: %d\n", nrow(AS)))

hr("1. PREMISE TEST — do two references agree on the same focal copy?")
cat("  Chance for a binary call with 2 references is 0.50.\n\n")
ag <- AS %>% group_by(focal, anchor, fg) %>% filter(n_distinct(ref) >= 2) %>%
  summarise(nref=n_distinct(ref),
            agree=n_distinct(rside)==1,
            majority=max(table(rside))/n(), .groups="drop")
print(as.data.frame(ag %>% group_by(focal) %>%
  summarise(copies=n(), med_refs=median(nref), frac_agree=round(mean(agree),3),
            med_majority=round(median(majority),3),
            p=signif(binom.test(sum(agree), n(), 0.5)$p.value,3), .groups="drop") %>%
  arrange(desc(frac_agree))), row.names=FALSE)
cat("\n  Well above 0.50 = the references define the same subgenomes and the\n")
cat("  shared-hybridisation premise is earned rather than assumed.\n")
cat("\n  by reference pair:\n")
print(as.data.frame(AS %>% group_by(focal, anchor, fg) %>%
  filter(n_distinct(ref)==2) %>%
  summarise(pair=paste(sort(sub("Drosera_","",ref)), collapse="+"),
            agree=n_distinct(rside)==1, .groups="drop") %>%
  group_by(pair) %>% summarise(n=n(), frac_agree=round(mean(agree),3),
                               .groups="drop") %>% arrange(desc(frac_agree))),
  row.names=FALSE)

hr("2. CONSTITUTION per focal species per reference")
res <- AS %>% group_by(focal, ref, anchor) %>%
  summarise(A=sum(rside=="A"), B=sum(rside=="B"), .groups="drop")
print(as.data.frame(res %>% group_by(focal, ref) %>%
  summarise(loci=n(), totA=sum(A), totB=sum(B),
            ratio=round(sum(A)/pmax(1,sum(B)),2), .groups="drop") %>%
  filter(loci>=30) %>% arrange(focal, desc(loci))), row.names=FALSE)
cat("\n  A focal species should give the SAME ratio against every reference.\n")

hr("3. consensus constitution — majority vote across references")
cons <- AS %>% group_by(focal, anchor, fg) %>%
  summarise(side=names(sort(table(rside), decreasing=TRUE))[1],
            nref=n_distinct(ref), unanimous=n_distinct(rside)==1, .groups="drop")
print(as.data.frame(cons %>% group_by(focal) %>%
  summarise(copies=n(), A=sum(side=="A"), B=sum(side=="B"),
            ratio=round(sum(side=="A")/pmax(1,sum(side=="B")),2), .groups="drop")),
  row.names=FALSE)
cat("\n  unanimous calls only:\n")
print(as.data.frame(cons %>% filter(unanimous | nref==1) %>% group_by(focal) %>%
  summarise(copies=n(), A=sum(side=="A"), B=sum(side=="B"),
            ratio=round(sum(side=="A")/pmax(1,sum(side=="B")),2), .groups="drop")),
  row.names=FALSE)

hr("4. per-locus composition, consensus calls")
loc <- cons %>% group_by(focal, anchor) %>%
  summarise(n=n(), A=sum(side=="A"), .groups="drop") %>% filter(n>=2)
for (f in CORE) {
  d <- loc %>% filter(focal==f)
  cat(sprintf("\n  %s — loci %d\n", sub("Drosera_","",f), nrow(d)))
  print(as.data.frame(d %>% group_by(n, A) %>% summarise(loci=n(), .groups="drop") %>%
    pivot_wider(names_from=A, values_from=loci, values_fill=0, names_prefix="A=")),
    row.names=FALSE)
  d3 <- d %>% filter(n==3)
  if (nrow(d3) >= 10) {
    p <- sum(d$A)/sum(d$n)
    cat(sprintf("    3-copy loci: %d | 2A:1B in %d (%.2f) | expected %.2f | p = %.3g\n",
      nrow(d3), sum(d3$A==2), mean(d3$A==2), 3*p^2*(1-p),
      binom.test(sum(d3$A==2), nrow(d3), 3*p^2*(1-p))$p.value))
  }
}
write_csv(AS, "trees/qc/crossref_fixed.csv")
write_csv(cons, "trees/qc/crossref_consensus.csv")
