#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"
CORE <- c("Drosera_binata","Drosera_paradoxa","Drosera_scorpioides")

k <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
     filter(!is.na(dS), dS>0, dS<5, codons>=100)
WI <- k %>% filter(sp1==sp2, sp1!=NEP) %>% transmute(anchor, sp=sp1, g1=seq1, g2=seq2, dS)
CR <- k %>% filter(sp1!=sp2, sp1!=NEP, sp2!=NEP) %>%
      mutate(a=pmin(sp1,sp2), b=pmax(sp1,sp2),
             ga=ifelse(sp1==a,seq1,seq2), gb=ifelse(sp1==a,seq2,seq1)) %>%
      select(anchor, a, b, ga, gb, dS)
cop <- bind_rows(CR %>% transmute(anchor, sp=a, g=ga), CR %>% transmute(anchor, sp=b, g=gb),
                 WI %>% transmute(anchor, sp, g=g1), WI %>% transmute(anchor, sp, g=g2)) %>%
       distinct()

hr("1. within-species dS — where the recent duplications sit")
cat("  The sweep must cross each species' recent peak without reaching its A/B split.\n\n")
print(as.data.frame(WI %>% group_by(sp) %>%
  summarise(pairs=n(), p10=round(quantile(dS,.10),3), p25=round(quantile(dS,.25),3),
            med=round(median(dS),3), p75=round(quantile(dS,.75),3), .groups="drop")),
  row.names=FALSE)

## single-linkage collapse of copies at threshold TT
clusters_at <- function(TT) {
  ed <- WI %>% filter(dS < TT)
  cop %>% group_by(anchor, sp) %>% group_modify(function(d, key) {
    ids <- d$g; lab <- setNames(seq_along(ids), ids)
    e <- ed[ed$anchor==key$anchor & ed$sp==key$sp, , drop=FALSE]
    e <- e[e$g1 %in% ids & e$g2 %in% ids, , drop=FALSE]
    if (nrow(e)) for (it in seq_along(ids)) for (j in seq_len(nrow(e))) {
      x <- lab[[e$g1[j]]]; y <- lab[[e$g2[j]]]
      if (x != y) { m <- min(x,y); lab[lab %in% c(x,y)] <- m }
    }
    tibble(g = ids, cl = unname(lab))
  }) %>% ungroup()
}

run_at <- function(TT) {
  cl <- clusters_at(TT)
  agg <- CR %>%
    inner_join(cl %>% rename(a=sp, ga=g, cla=cl), by=c("anchor","a","ga")) %>%
    inner_join(cl %>% rename(b=sp, gb=g, clb=cl), by=c("anchor","b","gb")) %>%
    group_by(anchor, a, b, cla, clb) %>%
    summarise(d = mean(dS), .groups="drop")
  g <- agg %>% group_by(anchor, a, b) %>%
    filter(n()==4, n_distinct(cla)==2, n_distinct(clb)==2) %>%
    mutate(ra=dense_rank(cla), rb=dense_rank(clb)) %>% ungroup() %>%
    select(anchor,a,b,ra,rb,d) %>% unite("cell", ra, rb) %>%
    pivot_wider(names_from=cell, values_from=d, names_prefix="m")
  need <- c("m1_1","m1_2","m2_1","m2_2")
  if (!all(need %in% names(g))) return(NULL)
  g %>% filter(if_all(all_of(need), ~!is.na(.))) %>%
    mutate(s1=m1_1+m2_2, s2=m1_2+m2_1,
           p1=ifelse(s1<=s2,m1_1,m1_2), p2=ifelse(s1<=s2,m2_2,m2_1),
           dd1=pmin(p1,p2), dd2=pmax(p1,p2), r=dd2/dd1, T=TT)
}

hr("2. threshold sweep (takes a couple of minutes)")
TS <- c(0.001, 0.15, 0.25, 0.35, 0.45, 0.60)
all <- bind_rows(lapply(TS, function(TT) {
  x <- run_at(TT); if (is.null(x)) return(NULL)
  x %>% group_by(a, b, T) %>% filter(n() >= 20) %>%
    summarise(loci=n(), med_d1=round(median(dd1),3), med_d2=round(median(dd2),3),
              ratio=round(median(r),3), .groups="drop")
}))
cat("  ratio here is the median of PER-LOCUS d2/d1 (the sec-13 bug is fixed).\n\n")
print(as.data.frame(all %>%
  mutate(pair=paste(sub("Dionaea_muscipula","DIO",sub("Drosera_","",a)),
                    sub("Dionaea_muscipula","DIO",sub("Drosera_","",b)), sep="/")) %>%
  select(pair, T, loci, med_d1, med_d2, ratio) %>%
  pivot_wider(names_from=T, values_from=c(loci, ratio), id_cols=pair)), row.names=FALSE)

hr("3. CONTROL — core trio must stay flat")
cat("  These have no recent duplications to merge. Drift = over-collapsing.\n\n")
print(as.data.frame(all %>% filter(a %in% CORE, b %in% CORE) %>%
  select(a,b,T,loci,ratio) %>% arrange(a,b,T)), row.names=FALSE)

hr("4. regia and capensis across the sweep")
print(as.data.frame(all %>%
  filter(grepl("regia|capensis", a) | grepl("regia|capensis", b)) %>%
  select(a,b,T,loci,med_d1,ratio) %>% arrange(a,b,T)), row.names=FALSE)
cat("\n  Ratio falling toward 1.27 as T rises = the elevation WAS recent duplicates,\n")
cat("  and that species shares the A/B substrate. Ratio flat and high = it does not.\n")

hr("5. Dio/regia d1 anomaly")
cat("  At T=0 this was 0.276, less than half every other Dionaea pair.\n\n")
print(as.data.frame(all %>% filter(a=="Dionaea_muscipula", b=="Drosera_regia") %>%
  select(T, loci, med_d1, med_d2, ratio)), row.names=FALSE)
cat("\n  d1 rising toward ~0.6 = it was regia's recent duplicates. Staying at 0.28\n")
cat("  = regia's donors really are far closer to Dionaea's than any other Drosera's.\n")
write_csv(all, "trees/qc/collapse_sweep.csv")
