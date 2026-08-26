#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); OFF <- 0.0230; BAND <- 0.08
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"; DIO <- "Dionaea_muscipula"
SPP <- c("Drosera_binata","Drosera_paradoxa","Drosera_scorpioides",
         "Drosera_regia","Drosera_capensis")
X <- c("chr1_sg1_s5","chr2_sg1_s3","chr3_sg1_s7","chr4_sg2_s12",
       "chr5_sg2_s16","chr6_sg2_s8","chr7_sg1_s1","chr8_sg1_s4")

k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS>0, dS<5, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
key <- tm %>% select(anchor=nep_gene, k1=tip, genome, chr) %>% distinct() %>%
       mutate(region=sub("-.*$","",anchor))
DD <- bind_rows(k %>% transmute(anchor, ga=seq1, gb=seq2, d=dS),
                k %>% transmute(anchor, ga=seq2, gb=seq1, d=dS))

hr("1. how separable is ortholog vs different-subgenome, per species pair?")
cat("  Assignment works only if cross-species ortholog dS sits well below\n")
cat("  within-species different-subgenome dS.\n\n")
sep <- bind_rows(lapply(SPP, function(s) {
  wi <- k$dS[k$sp1==s & k$sp2==s]
  bind_rows(lapply(setdiff(SPP, s), function(t){
    cr <- k$dS[(k$sp1==s & k$sp2==t) | (k$sp1==t & k$sp2==s)]
    if (length(cr) < 50) return(NULL)
    tibble(focal=sub("Drosera_","",s), ref=sub("Drosera_","",t),
           cross_med=round(median(cr),3), within_med=round(median(wi),3),
           gap=round(median(wi)/median(cr),2),
           sd_sep=round(abs(log(median(wi))-log(median(cr)))/0.435,2))
  }))
}))
print(as.data.frame(sep %>% arrange(desc(sd_sep))), row.names=FALSE)
cat("\n  sd_sep is the separation in log-SD units (0.435 from script 26).\n")
cat("  Above ~2 the per-copy assignment is reliable; below ~1.5 it is not.\n")

## per-copy delta for A/B labelling of the reference copies
dio <- key %>% filter(genome==DIO) %>% mutate(xy=ifelse(chr %in% X,"X","Y")) %>%
  select(anchor,k1,xy) %>% pivot_wider(names_from=xy, values_from=k1, values_fn=list) %>%
  filter(lengths(X)==1, lengths(Y)==1) %>%
  mutate(DX=unlist(X), DY=unlist(Y)) %>% select(anchor,DX,DY)
DEL <- key %>% inner_join(dio, by="anchor") %>%
  inner_join(DD %>% rename(DX=ga, DY=gb, dXY=d), by=c("anchor","DX","DY")) %>%
  left_join(DD %>% rename(k1=ga, DX=gb, dRX=d), by=c("anchor","k1","DX")) %>%
  left_join(DD %>% rename(k1=ga, DY=gb, dRY=d), by=c("anchor","k1","DY")) %>%
  filter(!is.na(dRX), !is.na(dRY), dXY>0.05) %>%
  transmute(anchor, k1, delta=(dRX-dRY)/dXY - OFF)

## classify each FOCAL copy by its nearest REF copy
assign_one <- function(FO, RE, minref=2) {
  f <- key %>% filter(genome==FO) %>% select(anchor, fg=k1, fchr=chr, region)
  r <- key %>% filter(genome==RE) %>% group_by(anchor) %>%
       filter(n_distinct(chr) >= minref) %>% ungroup() %>%
       select(anchor, rg=k1, rchr=chr)
  f %>% inner_join(r, by="anchor", relationship="many-to-many") %>%
    left_join(DD %>% rename(fg=ga, rg=gb, d=d), by=c("anchor","fg","rg")) %>%
    filter(!is.na(d)) %>%
    group_by(anchor, region, fg, fchr) %>%
    arrange(d, .by_group=TRUE) %>%
    summarise(best_rchr = rchr[1], d1 = d[1],
              d2 = if (n() > 1) d[2] else NA_real_,
              margin = if (n() > 1) (d[2]-d[1])/d[1] else NA_real_,
              .groups="drop") %>%
    mutate(focal=FO, ref=RE)
}

hr("2. per-copy assignment confidence")
cat("  margin = (2nd nearest - nearest)/nearest. Large = unambiguous.\n\n")
pairs <- expand.grid(f=SPP, r=SPP, stringsAsFactors=FALSE) %>% filter(f != r)
AS <- bind_rows(lapply(seq_len(nrow(pairs)), function(i)
        assign_one(pairs$f[i], pairs$r[i]))) %>% filter(!is.na(margin))
print(as.data.frame(AS %>% group_by(focal, ref) %>%
  summarise(copies=n(), loci=n_distinct(anchor), med_margin=round(median(margin),2),
            frac_clear=round(mean(margin>0.5),3), .groups="drop") %>%
  filter(copies>=100) %>% arrange(desc(frac_clear)) %>% head(12)), row.names=FALSE)

hr("3. CONSTITUTION per focal species, per reference — nothing pooled")
cat("  Each focal copy is labelled A or B by the delta of its nearest ref copy.\n")
cat("  Confident assignments only (margin > 0.5).\n\n")
LAB <- AS %>% filter(margin > 0.5) %>%
  left_join(key %>% filter(genome %in% SPP) %>% select(anchor, rg=k1, rgen=genome),
            by=c("anchor")) %>% filter(rgen == ref) %>% distinct() %>%
  left_join(DEL %>% rename(rg=k1, rdelta=delta), by=c("anchor","rg")) %>%
  filter(!is.na(rdelta)) %>%
  mutate(side = case_when(rdelta < -BAND ~ "A", rdelta > BAND ~ "B", TRUE ~ "amb")) %>%
  filter(side != "amb")
res <- LAB %>% group_by(focal, ref, anchor) %>%
  summarise(A=sum(side=="A"), B=sum(side=="B"), .groups="drop")
print(as.data.frame(res %>% group_by(focal, ref) %>%
  summarise(loci=n(), totA=sum(A), totB=sum(B),
            ratio=round(sum(A)/pmax(1,sum(B)),2), .groups="drop") %>%
  filter(loci >= 30) %>% arrange(focal, desc(loci))), row.names=FALSE)
cat("\n  Same focal species should give the same ratio against every reference.\n")
cat("  Divergent ratios = the reference frames disagree.\n")

hr("4. PREMISE TEST — do two references agree on the SAME copy?")
cat("  This is what pooling assumed and never checked.\n\n")
cmp <- LAB %>% select(focal, ref, anchor, fg, side) %>% distinct() %>%
  group_by(focal, anchor, fg) %>% filter(n_distinct(ref) >= 2) %>%
  summarise(nref=n_distinct(ref), agree=n_distinct(side)==1, .groups="drop")
print(as.data.frame(cmp %>% group_by(focal) %>%
  summarise(copies=n(), med_refs=median(nref),
            frac_agree=round(mean(agree),3), .groups="drop") %>%
  arrange(desc(frac_agree))), row.names=FALSE)
cat("\n  frac_agree near 1 = one shared set of subgenomes; pooling would have been\n")
cat("  safe and this per-species version is simply cleaner.\n")
cat("  frac_agree near 0.5 = references disagree; the premise fails and every\n")
cat("  constitution number above is frame-dependent.\n")
write_csv(LAB, "trees/qc/crossref_assignments.csv")
