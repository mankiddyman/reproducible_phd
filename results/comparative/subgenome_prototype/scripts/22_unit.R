#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE)
GS <- if (!is.na(commandArgs(TRUE)[1])) commandArgs(TRUE)[1] else
      file.path(dirname(BASE),"genespace","results","combBed.txt")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"; DIO <- "Dionaea_muscipula"; OFF <- 0.0230
X <- c("chr1_sg1_s5","chr2_sg1_s3","chr3_sg1_s7","chr4_sg2_s12",
       "chr5_sg2_s16","chr6_sg2_s8","chr7_sg1_s1","chr8_sg1_s4")

k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS>0, dS<5, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
cb <- read_tsv(GS, show_col_types=FALSE) %>% transmute(k1=id, genome, ord)
key <- tm %>% select(anchor=nep_gene, k1=tip, gene, genome, chr) %>% distinct() %>%
       mutate(region=sub("-.*$","",anchor)) %>%
       left_join(cb %>% rename(gene=k1), by=c("gene","genome"))
cat(sprintf("position join: %.4f\n", mean(!is.na(key$ord))))
DD <- bind_rows(k %>% transmute(anchor, ga=seq1, gb=seq2, d=dS),
                k %>% transmute(anchor, ga=seq2, gb=seq1, d=dS))
dio <- key %>% filter(genome==DIO) %>% mutate(xy=ifelse(chr %in% X,"X","Y")) %>%
  select(anchor,k1,xy) %>% pivot_wider(names_from=xy, values_from=k1, values_fn=list) %>%
  filter(lengths(X)==1, lengths(Y)==1) %>%
  mutate(DX=unlist(X), DY=unlist(Y)) %>% select(anchor,DX,DY)
D <- key %>% filter(genome != DIO, genome != NEP, !is.na(ord)) %>%
  inner_join(dio, by="anchor") %>%
  inner_join(DD %>% rename(DX=ga, DY=gb, dXY=d), by=c("anchor","DX","DY")) %>%
  left_join(DD %>% rename(k1=ga, DX=gb, dRX=d), by=c("anchor","k1","DX")) %>%
  left_join(DD %>% rename(k1=ga, DY=gb, dRY=d), by=c("anchor","k1","DY")) %>%
  filter(!is.na(dRX), !is.na(dRY), dXY>0.05) %>%
  mutate(delta=(dRX-dRY)/dXY - OFF)

hr("1. COLLISION TEST — is delta positionally structured WITHIN a segment?")
cat("  Collided segment: two delta plateaus at different positions along the chr.\n")
cat("  Real single-ancestry segment: flat.\n\n")
mr <- function(x,w=7){n<-length(x); if(n<w) return(rep(mean(x),n))
  vapply(seq_len(n),function(i)mean(x[max(1,i-w%/%2):min(n,i+w%/%2)]),numeric(1))}
sg <- D %>% group_by(genome, region, chr) %>% filter(n()>=20) %>%
  arrange(ord, .by_group=TRUE) %>%
  mutate(rl=mr(delta)) %>%
  summarise(copies=n(), med=round(median(delta),3),
            swing=round(max(rl)-min(rl),3),
            crosses0 = min(rl) < -0.05 & max(rl) > 0.05,
            ac1 = if (sd(delta)>0) round(cor(delta[-n()], delta[-1]),3) else 0,
            .groups="drop")
print(as.data.frame(sg %>% group_by(genome) %>%
  summarise(segments=n(), med_swing=round(median(swing),3),
            frac_crossing=round(mean(crosses0),3),
            med_ac1=round(median(ac1),3), .groups="drop") %>%
  arrange(desc(frac_crossing))), row.names=FALSE)
cat("\n  High frac_crossing = the unit is merging two ancestries. Predicted worst\n")
cat("  for paradoxa (6 chr) and scorpioides (4 chr); best for binata/regia (16-17).\n")

hr("2. paradoxa chr4_dom — the inverted region")
print(as.data.frame(sg %>% filter(genome=="Drosera_paradoxa") %>%
  mutate(region=sub("_dom","",region), chr=sub("_hap1","",chr)) %>%
  select(region, chr, copies, med, swing, crosses0) %>% arrange(region)),
  row.names=FALSE)
cat("\n  If chr4_dom segments swing across 0, it is a collision, not an inversion.\n")

hr("3. capensis — conversion or a clean recent WGD?")
cat("  Conversion is patchy: high block-to-block variance in own-copy dS.\n")
cat("  A clean WGD is uniform. Compare variance across species.\n\n")
own <- k %>% filter(sp1==sp2, sp1!=NEP) %>% transmute(anchor, sp=sp1, g1=seq1, dS) %>%
  left_join(key %>% select(anchor, g1=k1, chr, region), by=c("anchor","g1"))
bl <- own %>% filter(!is.na(chr)) %>% group_by(sp, region, chr) %>%
  filter(n()>=8) %>%
  summarise(pairs=n(), med=median(dS), .groups="drop")
print(as.data.frame(bl %>% group_by(sp) %>%
  summarise(blocks=n(), grand_med=round(median(med),3),
            block_sd=round(sd(med),3),
            cv=round(sd(med)/median(med),3),
            lo=round(min(med),3), hi=round(max(med),3), .groups="drop") %>%
  arrange(desc(cv))), row.names=FALSE)
cat("\n  capensis with the highest CV supports patchy conversion over a clean WGD.\n")
cat("\n  capensis blocks, sorted by own-copy dS:\n")
print(as.data.frame(bl %>% filter(sp=="Drosera_capensis") %>%
  mutate(region=sub("_dom","",region), chr=sub("_collapsed","",chr)) %>%
  select(region, chr, pairs, med) %>% arrange(med) %>% head(20)), row.names=FALSE)
cat("\n  A bimodal split — some blocks near 0.1, others near 0.7 — is conversion\n")
cat("  acting on part of the genome, with the rest retaining the ancient signal.\n")
write_csv(sg, "trees/qc/segment_positional.csv")
write_csv(bl, "trees/qc/block_owncopy_ds.csv")
