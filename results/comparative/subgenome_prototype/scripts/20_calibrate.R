#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); MINC <- 10
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"; DIO <- "Dionaea_muscipula"
X <- c("chr1_sg1_s5","chr2_sg1_s3","chr3_sg1_s7","chr4_sg2_s12",
       "chr5_sg2_s16","chr6_sg2_s8","chr7_sg1_s1","chr8_sg1_s4")
k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS>0, dS<5, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
key <- tm %>% select(anchor=nep_gene, k1=tip, genome, chr) %>% distinct() %>%
       mutate(region=sub("-.*$","",anchor))
DD <- bind_rows(k %>% transmute(anchor, ga=seq1, gb=seq2, d=dS),
                k %>% transmute(anchor, ga=seq2, gb=seq1, d=dS))
dio <- key %>% filter(genome==DIO) %>% mutate(xy=ifelse(chr %in% X,"X","Y")) %>%
       select(anchor,k1,xy) %>% pivot_wider(names_from=xy, values_from=k1, values_fn=list) %>%
       filter(lengths(X)==1, lengths(Y)==1) %>%
       mutate(DX=unlist(X), DY=unlist(Y)) %>% select(anchor,DX,DY)
D <- key %>% filter(genome != DIO) %>% inner_join(dio, by="anchor") %>%
  inner_join(DD %>% rename(DX=ga, DY=gb, dXY=d), by=c("anchor","DX","DY")) %>%
  left_join(DD %>% rename(k1=ga, DX=gb, dRX=d), by=c("anchor","k1","DX")) %>%
  left_join(DD %>% rename(k1=ga, DY=gb, dRY=d), by=c("anchor","k1","DY")) %>%
  filter(!is.na(dRX), !is.na(dRY), dXY>0.05) %>% mutate(delta=(dRX-dRY)/dXY)

hr("1. Nepenthes as the zero point")
np <- D %>% filter(genome==NEP)
cat(sprintf("  Nepenthes copies: %d\n", nrow(np)))
cat(sprintf("  median delta = %+.4f   [IQR %+.3f, %+.3f]\n",
            median(np$delta), quantile(np$delta,.25), quantile(np$delta,.75)))
cat(sprintf("  Wilcoxon vs 0: p = %.3g\n", wilcox.test(np$delta, mu=0)$p.value))
OFF <- median(np$delta)
cat(sprintf("\n  offset to subtract: %+.4f\n", OFF))
cat("  If this is near 0, the constitution ratios in 19 stand as printed.\n")
cat("  If clearly negative, every species was pushed toward A.\n")
cat("\n  by ancestral region:\n")
print(as.data.frame(np %>% group_by(region) %>%
  summarise(n=n(), med=round(median(delta),3), .groups="drop")), row.names=FALSE)

hr("2. constitution recomputed on the corrected axis")
seg <- D %>% filter(genome != NEP) %>% group_by(genome, region, chr) %>%
  summarise(copies=n(), med_delta=median(delta)-OFF, .groups="drop") %>%
  filter(copies >= MINC)
for (B in c(0.05, 0.08, 0.12)) {
  cc <- seg %>% mutate(side=case_when(med_delta < -B ~ "A",
                                      med_delta >  B ~ "B", TRUE ~ "amb"))
  cat(sprintf("\n  band +/- %.2f\n", B))
  print(as.data.frame(cc %>% group_by(genome) %>%
    summarise(A=sum(side=="A"), B=sum(side=="B"), amb=sum(side=="amb"),
              ratio=round(sum(side=="A")/pmax(1,sum(side=="B")),2), .groups="drop")),
    row.names=FALSE)
}
cat("\n  Ratios stable across bands = robust. Swinging = the band is doing the work.\n")
write_csv(seg, "trees/qc/segment_delta_calibrated.csv")
