#!/usr/bin/env Rscript
# regia_A is a mosaic: 743 loci where regia has 2 A copies and one is picked
# arbitrarily, with A1/A2 diverged at 0.268 (DR01, confirmed twice). A chimera
# of two divergent lineages is pulled basal. Rebuild using ONLY loci where
# regia has exactly ONE A copy, so no chimera is possible.
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
G  <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE)
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE)
lab <- SS %>% filter(label %in% c("A","B")) %>%
       select(genome, region, chr, mb_lo, mb_hi, label)
GL <- G %>% mutate(mb=mid/1e6) %>%
  inner_join(lab, by=c("genome","region","chr"), relationship="many-to-many") %>%
  filter(mb>=mb_lo, mb<=mb_hi) %>% group_by(tip) %>% slice(1) %>% ungroup()
ok <- GL %>% filter(genome=="Drosera_regia", label=="A") %>%
      count(locus) %>% filter(n==1) %>% pull(locus)
cat(sprintf("loci where regia has EXACTLY ONE A copy: %d\n", length(ok)))
LOC <- read_tsv("DR/locus_meta.tsv", show_col_types=FALSE) %>%
       group_by(tip) %>% slice(1) %>% ungroup()
P <- read_csv("fractionation_by_chrpair.csv", show_col_types=FALSE)
s1 <- P$retained_more == P$chrA
SIDE <- c(setNames(ifelse(s1,"A","B"), P$chrA), setNames(ifelse(s1,"B","A"), P$chrB))
dio <- LOC %>% filter(genome=="Dionaea_muscipula") %>%
  mutate(side=unname(SIDE[chr])) %>% filter(!is.na(side)) %>%
  transmute(locus, tip, tipname=paste0("Dionaea_", side))
nep <- LOC %>% filter(genome=="Nepenthes_gracilis") %>% group_by(locus) %>%
  slice(1) %>% ungroup() %>% transmute(locus, tip, tipname="Nepenthes")
readfa <- function(p) { x <- readLines(p, warn=FALSE); h <- grep("^>", x)
  nm <- sub("^>","",sub(" .*$","",x[h])); st <- h+1; en <- c(h[-1]-1, length(x))
  setNames(vapply(seq_along(h), function(i) paste0(x[st[i]:en[i]],collapse=""),
                  character(1)), nm) }
d <- GL %>% filter(locus %in% ok) %>%
  transmute(locus, tip, tipname=paste0(sub("Drosera_","",genome),"_",label)) %>%
  group_by(locus, tipname) %>% slice(1) %>% ungroup() %>%
  bind_rows(dio %>% filter(locus %in% ok), nep %>% filter(locus %in% ok))
KEEP <- c(paste0(c("regia","binata","paradoxa","scorpioides","capensis"),"_A"),
          "Dionaea_A","Nepenthes")
d <- d %>% filter(tipname %in% KEEP)
loci <- d %>% count(locus) %>% filter(n>=4) %>% pull(locus)
rows <- bind_rows(lapply(loci, function(l) {
  f <- file.path("DR/codon", paste0(l,".fna")); if (!file.exists(f)) return(NULL)
  s <- readfa(f); a <- d[d$locus==l,]; a <- a[a$tip %in% names(s),]
  if (nrow(a) < 4) return(NULL)
  L <- nchar(s[a$tip[1]]); if (L < 150 || L %% 3) return(NULL)
  tibble(locus=l, tipname=a$tipname, seq=unname(s[a$tip]), len=L) }))
Lx <- rows %>% distinct(locus, len) %>% arrange(locus)
mat <- setNames(vapply(KEEP, function(t) paste0(vapply(Lx$locus, function(l) {
  v <- rows$seq[rows$locus==l & rows$tipname==t]
  if (length(v)) v[1] else strrep("-", Lx$len[Lx$locus==l]) },
  character(1)), collapse=""), character(1)), KEEP)
writeLines(c(sprintf(" %d %d", length(mat), nchar(mat[1])),
             paste0(names(mat), "  ", mat)), "DR/tree/A_regia1copy.phy")
cat(sprintf("wrote DR/tree/A_regia1copy.phy: %d tips, %d loci, %d sites\n",
            length(mat), nrow(Lx), nchar(mat[1])))
cat("\nIF regia_A moves INSIDE Drosera here, the basal placement in A_full was\n")
cat("chimerism, not biology. IF it stays basal, the signal is real.\n")
