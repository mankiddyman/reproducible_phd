#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE)
GS <- file.path(dirname(BASE),"genespace","results","combBed.txt")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
DIO <- "Dionaea_muscipula"
PAIR <- c(chr1_sg1_s5="p1", chr1_sg2_s9="p1", chr2_sg1_s3="p2", chr2_sg2_s11="p2",
          chr3_sg1_s7="p3", chr3_sg2_s15="p3", chr4_sg1_s2="p4", chr4_sg2_s12="p4",
          chr5_sg1_s10="p5", chr5_sg2_s16="p5", chr6_sg1_s6="p6", chr6_sg2_s8="p6",
          chr7_sg1_s1="p7", chr7_sg2_s13="p7", chr8_sg1_s4="p8", chr8_sg2_s14="p8")

k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>% filter(!is.na(dS), dS>0, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
cb <- read_tsv(GS, show_col_types=FALSE) %>% transmute(gene=id, genome, ord, start, end)
key <- tm %>% select(anchor=nep_gene, k1=tip, gene, genome, chr) %>% distinct() %>%
       left_join(cb, by=c("gene","genome")) %>% mutate(pair=unname(PAIR[chr]))

D <- k %>% filter(sp1==DIO, sp2==DIO) %>% transmute(anchor, g1=seq1, g2=seq2, dS) %>%
  left_join(key %>% select(anchor, g1=k1, c1=chr, p1=pair, o1=ord, s1=start), by=c("anchor","g1")) %>%
  left_join(key %>% select(anchor, g2=k1, c2=chr, p2=pair, o2=ord, s2=start), by=c("anchor","g2")) %>%
  filter(!is.na(c1), !is.na(c2)) %>%
  mutate(cls = case_when(c1==c2 ~ "same chromosome",
                         p1==p2 ~ "same pair, different chr",
                         TRUE   ~ "DIFFERENT pair"))

hr("1. how are Dionaea's two copies arranged?")
print(as.data.frame(D %>% count(cls, name="pairs") %>%
  mutate(pct=round(100*pairs/sum(pairs),1))), row.names=FALSE)
cat("\n  v2 sec3.1 reported 16 same-chromosome and 191/906 cross-pair (21%).\n")
cat("  'DIFFERENT pair' means the two copies of one ancestral gene sit on\n")
cat("  chromosomes belonging to two different homeolog pairs — translocation.\n")

hr("2. dS by arrangement — is there a homeolog-depth subpopulation?")
cat("  Genuine homeologs sit near 0.525. Ancient dispersed paralogs sit much deeper.\n\n")
print(as.data.frame(D %>% group_by(cls) %>%
  summarise(n=n(), p05=round(quantile(dS,.05),3), p25=round(quantile(dS,.25),3),
            med=round(median(dS),3), p75=round(quantile(dS,.75),3),
            frac_near_homeolog=round(mean(dS>0.30 & dS<0.85),3), .groups="drop")),
  row.names=FALSE)
cat("\n  frac_near_homeolog = fraction inside the 0.30-0.85 homeolog window.\n")
cat("  A substantial fraction among same-chromosome pairs = translocated homeologs.\n")

hr("3. the candidates — same-chromosome pairs at homeolog depth")
cand <- D %>% filter(cls=="same chromosome", dS>0.30, dS<0.85) %>%
  mutate(gap_Mb = round(abs(s1-s2)/1e6,2), gap_genes = abs(o1-o2)) %>%
  select(anchor, chr=c1, dS, gap_Mb, gap_genes) %>% arrange(chr, gap_Mb)
cat(sprintf("  candidates: %d\n\n", nrow(cand)))
if (nrow(cand)) {
  print(as.data.frame(cand %>% mutate(dS=round(dS,3))), row.names=FALSE)
  cat("\n  by chromosome:\n")
  print(as.data.frame(cand %>% count(chr, name="n") %>% arrange(desc(n))), row.names=FALSE)
  cat("\n  Small gap_genes = tandem/local duplicate, NOT a translocated homeolog.\n")
  cat("  Large gap with homeolog-depth dS = the pattern you are describing.\n")
  cat(sprintf("\n  gap > 50 genes apart: %d of %d\n", sum(cand$gap_genes>50), nrow(cand)))
}

hr("4. do the candidates cluster positionally? (an HE tract would)")
if (nrow(cand) >= 5) {
  print(as.data.frame(cand %>% group_by(chr) %>% filter(n()>=2) %>%
    summarise(n=n(), span_Mb=round(diff(range(gap_Mb)),1),
              positions=paste(sort(round(gap_Mb,1)), collapse=", "), .groups="drop")),
    row.names=FALSE)
  cat("\n  Several candidates in one region = a genuine exchanged tract.\n")
  cat("  Scattered singletons = individual gene movements or noise.\n")
}

hr("5. DIFFERENT-pair copies — which chromosome pairs exchange?")
xp <- D %>% filter(cls=="DIFFERENT pair") %>% count(p1, p2, name="n") %>%
  filter(!is.na(p1), !is.na(p2)) %>% arrange(desc(n))
cat(sprintf("  distinct pair-combinations seen: %d\n\n", nrow(xp)))
print(as.data.frame(head(xp, 15)), row.names=FALSE)
cat("\n  Concentrated in a few combinations = real translocations between\n")
cat("  specific chromosomes. Spread evenly = orthology assignment noise.\n")

p <- D %>% ggplot(aes(dS)) + geom_histogram(bins=60, fill="#378ADD", colour="white") +
  geom_vline(xintercept=c(0.30,0.525,0.85), linetype=c("dashed","solid","dashed"),
             colour="firebrick") +
  facet_wrap(~cls, scales="free_y", ncol=1) + xlim(0,3) +
  labs(title="FIG53 - Dionaea copy pairs by chromosomal arrangement",
       subtitle="red line = 0.525, the homeolog depth; dashed = the 0.30-0.85 window",
       x="dS", y="pairs") + theme_minimal(10)
suppressWarnings(ggsave("trees/qc/FIG53_dionaea_arrangement.pdf", p, width=8, height=8))
write_csv(D %>% select(anchor, c1, c2, p1, p2, cls, dS), "trees/qc/dionaea_pair_arrangement.csv")
write_csv(cand, "trees/qc/dionaea_he_candidates.csv")
cat("\nWROTE: trees/qc/FIG53_dionaea_arrangement.pdf, dionaea_he_candidates.csv\n")
