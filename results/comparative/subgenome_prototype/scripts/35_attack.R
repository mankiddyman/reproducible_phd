#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); NPERM <- 200; BAND <- 0.08
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"; DIO <- "Dionaea_muscipula"
X <- c("chr1_sg1_s5","chr2_sg1_s3","chr3_sg1_s7","chr4_sg2_s12",
       "chr5_sg2_s16","chr6_sg2_s8","chr7_sg1_s1","chr8_sg1_s4")
spread <- function(v) round(max(v)/min(v), 3)

k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS>0, dS<3, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
key <- tm %>% select(anchor=nep_gene, k1=tip, genome, chr) %>% distinct() %>%
       mutate(region=sub("-.*$","",anchor))
RT <- read_csv("trees/qc/lineage_rates.csv", show_col_types=FALSE) %>%
      filter(dio_est=="dmin") %>% select(sp=species, rate)
DD <- bind_rows(k %>% transmute(anchor, ga=seq1, gb=seq2, d=dS),
                k %>% transmute(anchor, ga=seq2, gb=seq1, d=dS))

hr("ATTACK 1 — does the method collapse things that should NOT collapse?")
cat("  Same rates, same correction, applied to quantities that are known to\n")
cat("  differ. If everything collapses, the A/B result means nothing.\n\n")
seg <- read_csv("trees/qc/segment_delta_calibrated.csv", show_col_types=FALSE) %>%
  mutate(side=case_when(med_delta < -BAND ~ "A", med_delta > BAND ~ "B", TRUE ~ "amb")) %>%
  select(genome, region, chr, side)
P <- k %>% filter(sp1==sp2, sp1!=NEP) %>% transmute(anchor, sp=sp1, g1=seq1, g2=seq2, dS) %>%
  left_join(key %>% select(anchor, g1=k1, c1=chr, region), by=c("anchor","g1")) %>%
  left_join(key %>% select(anchor, g2=k1, c2=chr), by=c("anchor","g2")) %>%
  filter(!is.na(c1), !is.na(c2), c1!=c2) %>%
  left_join(seg %>% rename(sp=genome, c1=chr, s1=side), by=c("sp","region","c1")) %>%
  left_join(seg %>% rename(sp=genome, c2=chr, s2=side), by=c("sp","region","c2")) %>%
  filter(!is.na(s1), !is.na(s2), s1!="amb", s2!="amb") %>%
  mutate(cls=ifelse(s1==s2, paste0("same-",s1), "cross-AB")) %>%
  left_join(RT, by="sp")
ab   <- P %>% filter(cls=="cross-AB") %>% group_by(sp) %>%
        summarise(T=median(dS)/(2*rate[1]), .groups="drop")
sameA<- P %>% filter(cls=="same-A") %>% group_by(sp) %>%
        summarise(T=median(dS)/(2*rate[1]), .groups="drop")
sc <- tm %>% count(nep_gene, genome, name="kk") %>% filter(kk==1)
SPEC <- k %>% filter(sp1!=sp2, sp1!=NEP, sp2!=NEP, sp1!=DIO, sp2!=DIO) %>%
  semi_join(sc, by=c("anchor"="nep_gene","sp1"="genome")) %>%
  semi_join(sc, by=c("anchor"="nep_gene","sp2"="genome")) %>%
  mutate(a=pmin(sp1,sp2), b=pmax(sp1,sp2)) %>%
  left_join(RT %>% rename(a=sp, ra=rate), by="a") %>%
  left_join(RT %>% rename(b=sp, rb=rate), by="b") %>%
  group_by(a,b) %>% summarise(T=median(dS)/(ra[1]+rb[1]), .groups="drop")
cat(sprintf("  cross-AB depth   : %s   spread %.2fx\n",
    paste(round(sort(ab$T),3), collapse=" "), spread(ab$T)))
cat(sprintf("  same-A depth     : %s   spread %.2fx\n",
    paste(round(sort(sameA$T),3), collapse=" "), spread(sameA$T)))
cat(sprintf("  speciation depth : %s   spread %.2fx\n",
    paste(round(sort(SPEC$T),3), collapse=" "), spread(SPEC$T)))
cat("\n  Only cross-AB collapses. The other two use the SAME rates and the\n")
cat("  SAME correction and stay spread. The method does not force agreement.\n")

hr("ATTACK 2 — do the Dionaea X/Y labels earn their keep?")
cat("  Flip the X/Y assignment on random subsets of the 8 chromosome pairs,\n")
cat("  recompute segments and cross-AB, and see if the collapse survives.\n")
cat("  If a random orientation collapses just as well, the labels do nothing.\n\n")
dio <- key %>% filter(genome==DIO) %>% mutate(xy=ifelse(chr %in% X,"X","Y"),
                                              pair=sub("_sg.*","",chr)) %>%
  select(anchor, k1, xy, pair) %>%
  pivot_wider(names_from=xy, values_from=k1, values_fn=list) %>%
  filter(lengths(X)==1, lengths(Y)==1) %>%
  mutate(DX=unlist(X), DY=unlist(Y), pair=vapply(pair, function(z) z[1], character(1))) %>%
  select(anchor, DX, DY, pair)
BASEDEL <- key %>% filter(genome != DIO, genome != NEP) %>%
  inner_join(dio, by="anchor") %>%
  inner_join(DD %>% rename(DX=ga, DY=gb, dXY=d), by=c("anchor","DX","DY")) %>%
  left_join(DD %>% rename(k1=ga, DX=gb, dRX=d), by=c("anchor","k1","DX")) %>%
  left_join(DD %>% rename(k1=ga, DY=gb, dRY=d), by=c("anchor","k1","DY")) %>%
  filter(!is.na(dRX), !is.na(dRY), dXY>0.05) %>%
  transmute(anchor, k1, genome, chr, region, pair, delta=(dRX-dRY)/dXY - 0.0230)
PAIRS <- sort(unique(BASEDEL$pair))
run_orient <- function(flip) {
  d <- BASEDEL %>% mutate(delta = ifelse(pair %in% flip, -delta, delta))
  sg <- d %>% group_by(genome, region, chr) %>% filter(n()>=10) %>%
    summarise(md=median(delta), .groups="drop") %>%
    mutate(side=case_when(md < -BAND ~ "A", md > BAND ~ "B", TRUE ~ "amb")) %>%
    filter(side!="amb") %>% select(genome, region, chr, side)
  q <- k %>% filter(sp1==sp2, sp1!=NEP) %>% transmute(anchor, sp=sp1, g1=seq1, g2=seq2, dS) %>%
    left_join(key %>% select(anchor, g1=k1, c1=chr, region), by=c("anchor","g1")) %>%
    left_join(key %>% select(anchor, g2=k1, c2=chr), by=c("anchor","g2")) %>%
    filter(!is.na(c1), !is.na(c2), c1!=c2) %>%
    left_join(sg %>% rename(sp=genome, c1=chr, s1=side), by=c("sp","region","c1")) %>%
    left_join(sg %>% rename(sp=genome, c2=chr, s2=side), by=c("sp","region","c2")) %>%
    filter(!is.na(s1), !is.na(s2), s1!=s2) %>% left_join(RT, by="sp") %>%
    group_by(sp) %>% filter(n()>=25) %>%
    summarise(T=median(dS)/(2*rate[1]), .groups="drop")
  if (nrow(q) < 4) return(NA_real_)
  spread(q$T)
}
obs <- run_orient(character(0))
set.seed(1)
nul <- vapply(seq_len(NPERM), function(i)
  run_orient(sample(PAIRS, sample(1:7, 1))), numeric(1))
nul <- nul[!is.na(nul)]
cat(sprintf("  true orientation : spread %.3fx\n", obs))
cat(sprintf("  random flips     : median %.3fx   [%.3f - %.3f]   n = %d\n",
            median(nul), quantile(nul,.05), quantile(nul,.95), length(nul)))
cat(sprintf("  random orientations at least as tight: %.3f\n", mean(nul <= obs)))
cat("\n  A small p means the true X/Y orientation is special. A large p means\n")
cat("  any partition of Dionaea's chromosomes would have worked, and the A/B\n")
cat("  labels are not carrying information.\n")

hr("ATTACK 3 — how sensitive is the collapse to the fitted rates?")
cat("  Rate CIs from script 32: regia 0.210-0.395, Dionaea 0.285-0.420,\n")
cat("  paradoxa 1.059-1.240, scorpioides 1.077-1.253, capensis 0.912-1.176.\n\n")
raw <- P %>% filter(cls=="cross-AB") %>% group_by(sp) %>%
       summarise(raw=median(dS), n=n(), .groups="drop") %>% left_join(RT, by="sp")
CI <- tribble(~sp,~lo,~hi,
  "Drosera_regia",0.210,0.395, "Drosera_binata",1.000,1.000,
  "Drosera_paradoxa",1.059,1.240, "Drosera_scorpioides",1.077,1.253,
  "Drosera_capensis",0.912,1.176)
set.seed(2)
sw <- vapply(seq_len(2000), function(i) {
  r <- raw %>% left_join(CI, by="sp") %>%
       mutate(rr=runif(n(), lo, hi), T=raw/(2*rr))
  spread(r$T)
}, numeric(1))
cat(sprintf("  observed spread %.3fx\n", spread((raw$raw/(2*raw$rate)))))
cat(sprintf("  drawing rates uniformly from their CIs: median %.3fx  95%% %.3f-%.3f\n",
            median(sw), quantile(sw,.025), quantile(sw,.975)))
cat("\n  If the upper bound stays near 1.2x, the collapse is robust to rate\n")
cat("  uncertainty. If it reaches 1.6x+, the result is a property of one\n")
cat("  particular rate estimate.\n")

hr("ATTACK 4 — could we detect a real difference if there were one?")
cat("  Bootstrap each species' cross-AB median, then ask what spread we would\n")
cat("  see if the TRUE depths differed by 20% or 40%.\n\n")
pw <- P %>% filter(cls=="cross-AB")
bsd <- vapply(seq_len(1000), function(i) {
  q <- pw %>% group_by(sp) %>%
       summarise(T=median(sample(dS, n(), TRUE))/(2*rate[1]), .groups="drop")
  spread(q$T)
}, numeric(1))
cat(sprintf("  bootstrap spread under the observed data: median %.3f  95%% %.3f-%.3f\n",
            median(bsd), quantile(bsd,.025), quantile(bsd,.975)))
for (f in c(1.2, 1.4)) {
  sim <- vapply(seq_len(1000), function(i) {
    mult <- seq(1, f, length.out=nrow(raw))
    q <- pw %>% group_by(sp) %>%
         summarise(T=median(sample(dS, n(), TRUE))/(2*rate[1]), .groups="drop") %>%
         arrange(sp) %>% mutate(T=T*mult)
    spread(q$T)
  }, numeric(1))
  cat(sprintf("  if true depths spanned %.1fx: observed spread would be %.3f (median)\n",
              f, median(sim)))
}
cat("\n  If a 1.2x true difference would show as ~1.2x, the test has power and\n")
cat("  the observed 1.08x is meaningful. If a 1.2x difference would still read\n")
cat("  as ~1.1x, the collapse is not evidence of anything.\n")
