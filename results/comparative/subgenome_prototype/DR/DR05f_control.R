#!/usr/bin/env Rscript
# DR05f — is the regia+Dionaea signal the A/B LABELLING, or regia?
#
# The labels come from Dionaea: a segment is called A because it is closer to
# Dionaea_A. That biases d(focal_A, Dionaea_A) downward and favours exactly the
# answer DR05e returned.
#
# CONTROL: run the identical quartet with a CORE species in the focal slot.
#   (binata_A, Dionaea_A, paradoxa_A, Nepenthes)
# binata is unambiguously inside Drosera, so the true answer is focal+core.
#   binata+Dionaea also ~0.50  -> it is the labelling. DR05e collapses.
#   binata+Dionaea near 1/3    -> the bias is small and regia is genuinely
#                                 different. DR05e stands.
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
ALL <- c("regia","binata","paradoxa","scorpioides","capensis")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

LOC <- read_tsv("DR/locus_meta.tsv", show_col_types=FALSE) %>%
       group_by(tip) %>% slice(1) %>% ungroup()
K <- read_csv("DR/out/pairwise_ks.csv", show_col_types=FALSE) %>%
     filter(!is.na(dS), dS>=0, dS<5, codons>=100)
DD <- bind_rows(K %>% transmute(locus=anchor, a=seq1, b=seq2, d=dS),
                K %>% transmute(locus=anchor, a=seq2, b=seq1, d=dS))
dmap <- setNames(DD$d, paste(DD$locus, DD$a, DD$b))
gd <- function(l,x,y) ifelse(x==y, 0, unname(dmap[paste(l,x,y)]))
P <- read_csv("fractionation_by_chrpair.csv", show_col_types=FALSE)
s1 <- P$retained_more == P$chrA
SIDE <- c(setNames(ifelse(s1,"A","B"), P$chrA), setNames(ifelse(s1,"B","A"), P$chrB))
dio <- LOC %>% filter(genome=="Dionaea_muscipula") %>%
  mutate(side=unname(SIDE[chr])) %>% filter(!is.na(side)) %>%
  select(locus, tip, side) %>%
  pivot_wider(names_from=side, values_from=tip, values_fn=list) %>%
  filter(lengths(A)==1, lengths(B)==1) %>%
  transmute(locus, DioA=unlist(A), DioB=unlist(B))
nep <- LOC %>% filter(genome=="Nepenthes_gracilis") %>%
       group_by(locus) %>% slice(1) %>% ungroup() %>% transmute(locus, Nep=tip)
G  <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE)
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE)
lab <- SS %>% filter(label %in% c("A","B")) %>%
       select(genome, region, chr, mb_lo, mb_hi, label)
GL <- G %>% mutate(mb=mid/1e6) %>%
  inner_join(lab, by=c("genome","region","chr"), relationship="many-to-many") %>%
  filter(mb>=mb_lo, mb<=mb_hi) %>% group_by(tip) %>% slice(1) %>% ungroup() %>%
  mutate(sp=sub("Drosera_","",genome))
base <- dio %>% inner_join(nep, by="locus")

hr("EVERY species in the focal slot")
cat("  focal+Dionaea should be HIGH only for the species that really is\n")
cat("  sister to Dionaea. High for everyone = the labelling is doing it.\n\n")
res <- bind_rows(lapply(ALL, function(F) {
  oth <- setdiff(ALL, F)
  f <- GL %>% filter(sp==F) %>% select(locus, label, ftip=tip)
  o <- GL %>% filter(sp %in% oth) %>% select(locus, label, osp=sp, otip=tip)
  Q <- f %>% inner_join(o, by=c("locus","label"), relationship="many-to-many") %>%
       inner_join(base, by="locus") %>%
       mutate(Dio=ifelse(label=="A", DioA, DioB),
              S1=gd(locus,ftip,Dio)+gd(locus,otip,Nep),
              S2=gd(locus,ftip,otip)+gd(locus,Dio,Nep),
              S3=gd(locus,ftip,Nep)+gd(locus,Dio,otip)) %>%
       filter(!is.na(S1), !is.na(S2), !is.na(S3))
  if (!nrow(Q)) return(NULL)
  b <- apply(Q[,c("S1","S2","S3")],1,which.min)
  sr <- t(apply(Q[,c("S1","S2","S3")],1,sort))
  Q <- Q %>% mutate(vote=c("focal+Dionaea","focal+other","neither")[b],
                    gap=(sr[,2]-sr[,1])/rowMeans(cbind(S1,S2,S3))) %>%
       filter(gap > 0.02)
  L <- Q %>% count(locus, vote) %>% group_by(locus) %>%
    slice_max(n, n=1, with_ties=TRUE) %>%
    summarise(vote=if (n()==1) vote else NA_character_, .groups="drop") %>%
    filter(!is.na(vote))
  tb <- table(factor(L$vote, levels=c("focal+Dionaea","focal+other","neither")))
  bt <- binom.test(tb[["focal+Dionaea"]], sum(tb), 1/3)
  tibble(focal=F, loci=sum(tb),
         focal_Dionaea=round(tb[["focal+Dionaea"]]/sum(tb),3),
         focal_other=round(tb[["focal+other"]]/sum(tb),3),
         neither=round(tb[["neither"]]/sum(tb),3),
         CI=sprintf("[%.3f,%.3f]", bt$conf.int[1], bt$conf.int[2]),
         p=signif(bt$p.value,3))
}))
print(as.data.frame(res %>% arrange(desc(focal_Dionaea))), row.names=FALSE)
write_csv(res, "DR/out/DR05f_control.csv")

hr("VERDICT")
rg <- res$focal_Dionaea[res$focal=="regia"]
ot <- res$focal_Dionaea[res$focal!="regia"]
cat(sprintf("  regia   : %.3f\n  others  : %s\n  gap     : %+.3f\n",
            rg, paste(sprintf("%.3f", ot), collapse=", "), rg - max(ot)))
if (rg > max(ot) + 0.08) {
  cat("\n  regia is CLEARLY HIGHEST. The labelling bias affects all species\n")
  cat("  equally, so the excess for regia is real. DR05e STANDS.\n")
} else if (rg <= max(ot) + 0.03) {
  cat("\n  regia is NOT distinguishable from species known to be inside\n")
  cat("  Drosera. The signal is the LABELLING, not regia. DR05e COLLAPSES.\n")
} else {
  cat("\n  regia is somewhat higher but not decisively. Report as suggestive.\n")
}
cat("\n  Any elevation shared by ALL species is the label bias, and its size\n")
cat("  is measured by the non-regia rows.\n")
