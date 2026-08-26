#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); NPERM <- 499; OFF <- 0.0230
GS <- file.path(dirname(BASE),"genespace","results","combBed.txt")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"; DIO <- "Dionaea_muscipula"
X <- c("chr1_sg1_s5","chr2_sg1_s3","chr3_sg1_s7","chr4_sg2_s12",
       "chr5_sg2_s16","chr6_sg2_s8","chr7_sg1_s1","chr8_sg1_s4")

k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS>0, dS<5, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
cb <- read_tsv(GS, show_col_types=FALSE) %>% transmute(gene=id, genome, ord)
key <- tm %>% select(anchor=nep_gene, k1=tip, gene, genome, chr) %>% distinct() %>%
       mutate(region=sub("-.*$","",anchor)) %>% left_join(cb, by=c("gene","genome"))
DD <- bind_rows(k %>% transmute(anchor, ga=seq1, gb=seq2, d=dS),
                k %>% transmute(anchor, ga=seq2, gb=seq1, d=dS))
dio <- key %>% filter(genome==DIO) %>% mutate(xy=ifelse(chr %in% X,"X","Y")) %>%
  select(anchor,k1,xy) %>% pivot_wider(names_from=xy, values_from=k1, values_fn=list) %>%
  filter(lengths(X)==1, lengths(Y)==1) %>%
  mutate(DX=unlist(X), DY=unlist(Y)) %>% select(anchor,DX,DY)
D <- key %>% filter(genome != DIO, genome != NEP) %>% inner_join(dio, by="anchor") %>%
  inner_join(DD %>% rename(DX=ga, DY=gb, dXY=d), by=c("anchor","DX","DY")) %>%
  left_join(DD %>% rename(k1=ga, DX=gb, dRX=d), by=c("anchor","k1","DX")) %>%
  left_join(DD %>% rename(k1=ga, DY=gb, dRY=d), by=c("anchor","k1","DY")) %>%
  filter(!is.na(dRX), !is.na(dRY), dXY>0.05) %>% mutate(delta=(dRX-dRY)/dXY - OFF)

hr("1. positional structure, WITH a permutation null")
mr <- function(x,w=7){n<-length(x); if(n<w) return(rep(mean(x),n))
  vapply(seq_len(n),function(i)mean(x[max(1,i-w%/%2):min(n,i+w%/%2)]),numeric(1))}
set.seed(1)
pos <- D %>% filter(!is.na(ord)) %>% group_by(genome, region, chr) %>%
  filter(n() >= 25) %>% arrange(ord, .by_group=TRUE) %>%
  group_modify(function(d,key){
    v <- d$delta; o_sw <- diff(range(mr(v))); o_ac <- cor(v[-length(v)], v[-1])
    n_sw <- vapply(seq_len(NPERM), function(i) diff(range(mr(sample(v)))), numeric(1))
    n_ac <- vapply(seq_len(NPERM), function(i){s<-sample(v); cor(s[-length(s)], s[-1])}, numeric(1))
    tibble(n=length(v), swing=o_sw, p_swing=(1+sum(n_sw>=o_sw))/(NPERM+1),
           ac1=o_ac, p_ac=(1+sum(n_ac>=o_ac))/(NPERM+1))
  }) %>% ungroup()
print(as.data.frame(pos %>% group_by(genome) %>%
  summarise(segs=n(), med_swing=round(median(swing),3),
            null_beaten=round(mean(p_swing<0.05),3),
            ac_sig=round(mean(p_ac<0.05),3), .groups="drop")), row.names=FALSE)
cat("\n  null_beaten near 0 = the swing is noise and collision is undetectable here.\n")
cat("  Meaningfully above 0.05 = real positional structure within segments.\n")

hr("2. MINC sweep — is the constitution threshold-driven?")
for (M in c(5, 10, 20, 30)) {
  cc <- D %>% group_by(genome, region, chr) %>% filter(n() >= M) %>%
    summarise(med=median(delta), .groups="drop") %>%
    mutate(side=case_when(med < -0.08 ~ "A", med > 0.08 ~ "B", TRUE ~ "amb"))
  cat(sprintf("\n  MINC = %d\n", M))
  print(as.data.frame(cc %>% group_by(genome) %>%
    summarise(segs=n(), A=sum(side=="A"), B=sum(side=="B"),
              ratio=round(sum(side=="A")/pmax(1,sum(side=="B")),2), .groups="drop")),
    row.names=FALSE)
}
cat("\n  Ratio stable across MINC = real. Climbing as MINC rises = small B segments\n")
cat("  were being dropped, exactly the fractionation bias from v2 5.\n")

hr("3. CONVERSION TEST — within-species pairs by same vs opposite side")
cat("  Same-side pair = duplicate within one ancestral subgenome -> recent.\n")
cat("  Cross-side pair = A vs B homeolog -> ancient.\n")
cat("  Conversion = cross-side pairs pulled DOWN toward the same-side mode.\n\n")
seg <- D %>% group_by(genome, region, chr) %>% filter(n() >= 10) %>%
  summarise(med=median(delta), .groups="drop") %>%
  mutate(side=case_when(med < -0.08 ~ "A", med > 0.08 ~ "B", TRUE ~ "amb")) %>%
  filter(side != "amb")
wi <- k %>% filter(sp1==sp2, sp1!=NEP) %>% transmute(anchor, genome=sp1, g1=seq1, g2=seq2, dS) %>%
  left_join(key %>% select(anchor, g1=k1, genome, c1=chr, region), by=c("anchor","g1","genome")) %>%
  left_join(key %>% select(anchor, g2=k1, genome, c2=chr), by=c("anchor","g2","genome")) %>%
  left_join(seg %>% select(genome, region, c1=chr, s1=side), by=c("genome","region","c1")) %>%
  left_join(seg %>% select(genome, region, c2=chr, s2=side), by=c("genome","region","c2")) %>%
  filter(!is.na(s1), !is.na(s2)) %>%
  mutate(cls = ifelse(s1==s2, paste0("same-", s1), "cross"))
print(as.data.frame(wi %>% group_by(genome, cls) %>%
  summarise(pairs=n(), p10=round(quantile(dS,.10),3), med=round(median(dS),3),
            p90=round(quantile(dS,.90),3), .groups="drop") %>%
  arrange(genome, cls)), row.names=FALSE)
cat("\n  cross-side pairs with dS below the same-side median = converted homeologs:\n")
print(as.data.frame(wi %>% group_by(genome) %>%
  summarise(cross=sum(cls=="cross"),
            same_med=median(dS[cls!="cross"]),
            frac_cross_low=round(mean(dS[cls=="cross"] < median(dS[cls!="cross"])),3),
            .groups="drop")), row.names=FALSE)
write_csv(wi %>% select(anchor, genome, region, c1, c2, s1, s2, cls, dS),
          "trees/qc/pair_side_classes.csv")
