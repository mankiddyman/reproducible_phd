#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); NPERM <- 999; BAND <- 0.08
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
CORE <- c("Drosera_binata","Drosera_paradoxa","Drosera_scorpioides")
REFS <- c(CORE, "Drosera_regia")

tm  <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
key <- tm %>% select(anchor=nep_gene, k1=tip, genome, chr) %>% distinct() %>%
       mutate(region=sub("-.*$","",anchor))
seg <- read_csv("trees/qc/segment_delta_calibrated.csv", show_col_types=FALSE) %>%
  mutate(side=case_when(med_delta < -BAND ~ "A", med_delta > BAND ~ "B", TRUE ~ "amb")) %>%
  select(genome, region, chr, side)
AS <- read_csv("trees/qc/crossref_fixed.csv", show_col_types=FALSE)
cat("crossref_fixed columns:", paste(names(AS), collapse=", "), "\n")

hr("1. do the two focal copies ever land on the SAME reference copy?")
cat("  If nearest-neighbour matching ALWAYS gives distinct partners, the\n")
cat("  anti-correlation is structural, not biological.\n\n")
two <- AS %>% group_by(focal, ref, anchor) %>% filter(n()==2) %>%
  summarise(same_ref_chr = n_distinct(best_rchr)==1,
            same_side    = n_distinct(rside)==1, .groups="drop")
print(as.data.frame(two %>% group_by(focal) %>%
  summarise(loci=n(), frac_same_ref=round(mean(same_ref_chr),3),
            frac_same_side=round(mean(same_side),3), .groups="drop")), row.names=FALSE)
cat("\n  frac_same_ref near 0 = the matching forces distinct partners.\n")

hr("2. PERMUTATION — assign each focal copy a RANDOM reference copy")
cat("  Same loci, same reference copies available, same segment labels.\n")
cat("  Only the nearest-neighbour choice is replaced by a random one.\n")
cat("  If A=0 stays depleted, the depletion is mechanical.\n\n")
refpool <- key %>% filter(genome %in% REFS) %>%
  left_join(seg %>% rename(genome=genome, chr=chr), by=c("genome","region","chr")) %>%
  filter(!is.na(side), side!="amb") %>% select(anchor, ref=genome, rchr=chr, side)
foc <- AS %>% select(focal, ref, anchor, fg) %>% distinct()
obs <- AS %>% group_by(focal, anchor, fg) %>%
  summarise(side=names(sort(table(rside), decreasing=TRUE))[1], .groups="drop") %>%
  group_by(focal, anchor) %>% filter(n()==2) %>%
  summarise(A=sum(side=="A"), .groups="drop")
o <- obs %>% group_by(focal) %>%
  summarise(loci=n(), A0=sum(A==0), A1=sum(A==1), A2=sum(A==2),
            pA=(sum(A))/(2*n()), .groups="drop") %>%
  mutate(A0_indep=round((1-pA)^2*loci,1))
print(as.data.frame(o), row.names=FALSE)

set.seed(1)
nul <- vapply(seq_len(NPERM), function(i) {
  r <- foc %>% inner_join(refpool, by=c("anchor","ref"),
                          relationship="many-to-many") %>%
    group_by(focal, anchor, fg) %>% slice_sample(n=1) %>%
    group_by(focal, anchor, fg) %>%
    summarise(side=side[1], .groups="drop") %>%
    group_by(focal, anchor) %>% filter(n()==2) %>%
    summarise(A=sum(side=="A"), .groups="drop") %>%
    group_by(focal) %>% summarise(f0=mean(A==0), .groups="drop")
  setNames(r$f0, r$focal)[CORE]
}, numeric(length(CORE)))
cat("\n  fraction of 2-copy loci with A=0:\n")
print(as.data.frame(tibble(focal=CORE,
  observed = round(vapply(CORE, function(f){d<-obs%>%filter(focal==f); mean(d$A==0)}, numeric(1)),4),
  null_median = round(apply(nul,1,median,na.rm=TRUE),4),
  null_lo = round(apply(nul,1,quantile,.025,na.rm=TRUE),4),
  null_hi = round(apply(nul,1,quantile,.975,na.rm=TRUE),4))), row.names=FALSE)
cat("\n  Observed BELOW the random-assignment null = real signal.\n")
cat("  Observed INSIDE the null = the A=0 depletion in script 29 is an\n")
cat("  artefact of nearest-neighbour matching and must be withdrawn.\n")
