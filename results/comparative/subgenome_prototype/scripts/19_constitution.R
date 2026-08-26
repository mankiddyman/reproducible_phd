#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"; DIO <- "Dionaea_muscipula"; MINC <- 10; BAND <- 0.08
X <- c("chr1_sg1_s5","chr2_sg1_s3","chr3_sg1_s7","chr4_sg2_s12",
       "chr5_sg2_s16","chr6_sg2_s8","chr7_sg1_s1","chr8_sg1_s4")
PLOIDY <- c(Dionaea_muscipula=2, Drosera_regia=3, Drosera_binata=3,
            Drosera_paradoxa=3, Drosera_scorpioides=3, Drosera_capensis=6)

k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS>0, dS<5, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
key <- tm %>% select(anchor=nep_gene, k1=tip, genome, chr) %>% distinct() %>%
       mutate(region = sub("-.*$","",anchor))
DD <- bind_rows(k %>% transmute(anchor, ga=seq1, gb=seq2, d=dS),
                k %>% transmute(anchor, ga=seq2, gb=seq1, d=dS))

## delta axis from Dionaea's homeolog pair at each locus
dio <- key %>% filter(genome==DIO) %>% mutate(xy=ifelse(chr %in% X,"X","Y")) %>%
       select(anchor, k1, xy) %>%
       pivot_wider(names_from=xy, values_from=k1, values_fn=list) %>%
       filter(lengths(X)==1, lengths(Y)==1) %>%
       mutate(DX=unlist(X), DY=unlist(Y)) %>% select(anchor, DX, DY)
D <- key %>% filter(genome != DIO, genome != NEP) %>%
  inner_join(dio, by="anchor") %>%
  inner_join(DD %>% rename(DX=ga, DY=gb, dXY=d), by=c("anchor","DX","DY")) %>%
  left_join(DD %>% rename(k1=ga, DX=gb, dRX=d), by=c("anchor","k1","DX")) %>%
  left_join(DD %>% rename(k1=ga, DY=gb, dRY=d), by=c("anchor","k1","DY")) %>%
  filter(!is.na(dRX), !is.na(dRY), dXY > 0.05) %>%
  mutate(delta = (dRX - dRY)/dXY)

hr("1. segment inventory — does segment count per region match ploidy?")
cat("  segment = species x ancestral region x chromosome. This is the unit\n")
cat("  that survives mosaicism; a chromosome spanning 2 regions gives 2 segments.\n\n")
seg <- D %>% group_by(genome, region, chr) %>%
  summarise(copies=n(), med_delta=median(delta), .groups="drop")
print(as.data.frame(seg %>% filter(copies >= MINC) %>%
  group_by(genome, region) %>% summarise(segments=n(), .groups="drop") %>%
  group_by(genome) %>%
  summarise(regions=n(), med_segments=median(segments),
            rng=paste0(min(segments),"-",max(segments)), .groups="drop") %>%
  mutate(expected = PLOIDY[genome])), row.names=FALSE)
cat("\n  Compare to v2 8.3 major-tract counts: binata 3, regia 3, paradoxa 3,\n")
cat("  scorpioides 3, capensis 5. Shortfall = fractionation or the MINC filter.\n")

hr("2. mosaic chromosomes — same chromosome, different region, different side")
mos <- seg %>% filter(copies >= MINC) %>% group_by(genome, chr) %>%
  filter(n_distinct(region) >= 2) %>%
  summarise(regions=n(), span=round(max(med_delta)-min(med_delta),3),
            detail=paste(sprintf("%s:%+.2f", sub("_dom","",region), med_delta),
                         collapse=" "), .groups="drop") %>%
  filter(span > 0.15) %>% arrange(desc(span))
cat(sprintf("  chromosomes in >=2 regions with delta span > 0.15: %d\n\n", nrow(mos)))
print(as.data.frame(mos %>% mutate(genome=sub("Drosera_","D_",genome))), row.names=FALSE)
cat("\n  These are exactly the segments that chromosome-level pooling averages out.\n")

hr("3. CONSTITUTION — A:B segment counts per species per region")
cc <- seg %>% filter(copies >= MINC) %>%
  mutate(side = case_when(med_delta < -BAND ~ "A", med_delta > BAND ~ "B", TRUE ~ "amb"))
print(as.data.frame(cc %>% group_by(genome, region) %>%
  summarise(A=sum(side=="A"), B=sum(side=="B"), amb=sum(side=="amb"),
            .groups="drop") %>%
  mutate(call=paste0(A,":",B, ifelse(amb>0, paste0(" (+",amb,"?)"), ""))) %>%
  select(genome, region, call) %>%
  pivot_wider(names_from=region, values_from=call)), row.names=FALSE)

hr("4. modal constitution per species")
print(as.data.frame(cc %>% group_by(genome, region) %>%
  summarise(A=sum(side=="A"), B=sum(side=="B"), .groups="drop") %>%
  group_by(genome) %>%
  summarise(regions=n(), tot_A=sum(A), tot_B=sum(B),
            ratio=round(sum(A)/pmax(1,sum(B)),2),
            modal=names(sort(table(paste0(A,":",B)), decreasing=TRUE))[1],
            .groups="drop")), row.names=FALSE)
cat("\n  Counting SEGMENTS not copies: a segment is present or absent, so this is\n")
cat("  not skewed by the 34-41% gene-retention asymmetry (v2 5).\n")
cat("  regia should return 2:1 if 18 was right at the chromosome level too.\n")

hr("5. regia — segment level vs the chromosome-level version in 18")
print(as.data.frame(seg %>% filter(genome=="Drosera_regia", copies >= MINC) %>%
  mutate(chr=sub("_collapsed","",chr), region=sub("_dom","",region)) %>%
  select(region, chr, copies, med_delta) %>% arrange(region, med_delta)),
  row.names=FALSE)

p <- seg %>% filter(copies >= MINC) %>% mutate(genome=sub("Drosera_","D_",genome)) %>%
  ggplot(aes(med_delta)) + geom_histogram(binwidth=0.05, fill="#7F77DD", colour="white") +
  geom_vline(xintercept=c(-BAND,BAND), linetype="dashed", colour="firebrick") +
  facet_wrap(~genome, scales="free_y") +
  labs(title="FIG51 - segment placement on the Dionaea A/B axis",
       subtitle="one point per species x region x chromosome; bimodal = two ancestries",
       x="median delta", y="segments") + theme_minimal(10)
ggsave("trees/qc/FIG51_segments.pdf", p, width=10, height=6)
write_csv(seg, "trees/qc/segment_delta.csv")
cat("\nWROTE: trees/qc/segment_delta.csv  FIG51_segments.pdf\n")
