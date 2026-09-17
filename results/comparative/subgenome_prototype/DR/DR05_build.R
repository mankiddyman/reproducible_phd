#!/usr/bin/env Rscript
# ============================================================================
# DR05a — BUILD THE PER-SUBGENOME SUPERMATRICES
#
# THE QUESTION DR03b COULD NOT ANSWER
#   Does regia branch inside Drosera, or outside it? DR03b's two label sources
#   disagreed (-0.172 vs +0.108) and a rate change in regia explains the
#   segment-label result without any topological cause. Distance ratios cannot
#   separate those. TOPOLOGY can: IQ-TREE fits branch lengths as free
#   parameters, so lineage rate variation IS the branch lengths and does not
#   distort the tree shape.
#
# 13 TIPS
#   5 Drosera x {A, B} + Dionaea_A + Dionaea_B + Nepenthes
#   NOT 15: DR03b tested whether A1 and A2 are alignable across species and
#   found they are not (the rAA bimodality that motivated it turned out to be
#   label noise -- rBB was equally bimodal, dBIC 114.6 vs 96.6).
#
#   CHIMERISM, stated plainly: where a species has two A copies at a locus we
#   take one, so its A sequence is a mosaic of A1 and A2 across loci. Under
#   INDEPENDENT A-duplications that is harmless -- both copies coalesce with
#   other species at the same node. It would matter if the duplication were
#   shared, which DR03b argues against.
#
# THE MISLABEL PROBLEM, AND WHY NO THRESHOLD IS IMPOSED
#   DR03's mixture fit put 6-17% of loci in a component at ratio ~1. That is
#   an AGGREGATE, not a list -- it cannot say WHICH loci. And the component is
#   not necessarily error:
#     - mislabelled copy      -> would be error
#     - incomplete lineage sorting -> BIOLOGY. The A/B split (T~0.57) and the
#       regia speciation (T~0.44) are close, which is exactly the high-ILS
#       regime. Filtering these out would REMOVE the loci that disagree with
#       the majority and bias the tree.
#     - gene conversion between a species' own A and B copies
#
#   So instead of a cut: BOOTSTRAP each segment's own genes to get a per-
#   segment flip probability, then report the tree at full data AND at high-
#   stability segments. If the topology is identical, the mislabel rate is
#   irrelevant to the conclusion. If it changes, that IS the finding.
#
# IN   DR/out/DR02_gene_delta.csv, DR/out/DR02_segments.csv,
#      DR/locus_meta.tsv, DR/codon/*.fna
# OUT  DR/tree/A_full.phy, B_full.phy, A_stable.phy, B_stable.phy
#      DR/tree/genes/<locus>_{A,B}.fna   for the ASTRAL gene trees
#      DR/out/DR05_segment_stability.csv, DR05_matrix_stats.csv
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
set.seed(1); NB <- 1000L; MINGENE <- 8L
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
flush.console()

hr("1. SEGMENT STABILITY — bootstrap, not a threshold")
cat("  Resample each segment's genes with replacement, recompute the median\n")
cat("  delta, count how often the A/B call flips. That is a per-segment error\n")
cat("  rate derived FROM THE DATA, not a cut chosen by hand.\n\n")
G  <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE)
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE)
lab <- SS %>% filter(label %in% c("A","B")) %>%
       select(genome, region, chr, mb_lo, mb_hi, label, segment)
GL <- G %>% mutate(mb=mid/1e6) %>%
  inner_join(lab, by=c("genome","region","chr"), relationship="many-to-many") %>%
  filter(mb>=mb_lo, mb<=mb_hi) %>% group_by(tip) %>% slice(1) %>% ungroup()

STAB <- GL %>% group_by(segment, genome, label) %>%
  summarise(n=n(), med=median(delta), .groups="drop_last") %>% ungroup() %>%
  left_join(GL %>% group_by(segment) %>%
    summarise(flip = { d <- delta
      if (length(d) < 3) NA_real_ else {
        b <- replicate(NB, median(sample(d, replace=TRUE)))
        obs <- median(d); mean(sign(b) != sign(obs)) } }, .groups="drop"),
    by="segment")
cat("  flip probability across segments:\n")
print(summary(STAB$flip))
cat(sprintf("\n  segments with flip < 0.05: %d of %d (%.0f%%)\n",
            sum(STAB$flip < 0.05, na.rm=TRUE), nrow(STAB),
            100*mean(STAB$flip < 0.05, na.rm=TRUE)))
print(as.data.frame(STAB %>% group_by(genome) %>%
  summarise(segments=n(), median_flip=round(median(flip, na.rm=TRUE),3),
            stable=sum(flip<0.05, na.rm=TRUE), .groups="drop")), row.names=FALSE)
write_csv(STAB, "DR/out/DR05_segment_stability.csv")
GL <- GL %>% left_join(STAB %>% select(segment, flip), by="segment") %>%
      mutate(stable = !is.na(flip) & flip < 0.05)
cat(sprintf("\n  genes on stable segments: %d of %d (%.0f%%)\n",
            sum(GL$stable), nrow(GL), 100*mean(GL$stable)))
flush.console()

hr("2. ASSEMBLE per-locus, per-subgenome tip sequences")
LOC <- read_tsv("DR/locus_meta.tsv", show_col_types=FALSE) %>%
       group_by(tip) %>% slice(1) %>% ungroup()
P <- read_csv("fractionation_by_chrpair.csv", show_col_types=FALSE)
sg1_dom <- P$retained_more == P$chrA
SIDE <- c(setNames(ifelse(sg1_dom,"A","B"), P$chrA),
          setNames(ifelse(sg1_dom,"B","A"), P$chrB))
dio <- LOC %>% filter(genome=="Dionaea_muscipula") %>%
  mutate(side=unname(SIDE[chr])) %>% filter(!is.na(side)) %>%
  transmute(locus, tip, tipname=paste0("Dionaea_", side))
nep <- LOC %>% filter(genome=="Nepenthes_gracilis") %>%
  group_by(locus) %>% slice(1) %>% ungroup() %>%
  transmute(locus, tip, tipname="Nepenthes")
readfa <- function(p) {
  x <- readLines(p, warn=FALSE); h <- grep("^>", x)
  nm <- sub("^>", "", sub(" .*$", "", x[h]))
  st <- h+1; en <- c(h[-1]-1, length(x))
  setNames(vapply(seq_along(h), function(i)
    paste0(x[st[i]:en[i]], collapse=""), character(1)), nm)
}
TIPS <- c(paste0(rep(sub("Drosera_","",DROS), each=2), "_", c("A","B")),
          "Dionaea_A","Dionaea_B","Nepenthes")
cat(sprintf("  %d tips: %s\n", length(TIPS), paste(TIPS, collapse=", ")))

build <- function(use_stable) {
  d <- if (use_stable) GL %>% filter(stable) else GL
  d <- d %>% transmute(locus, tip,
        tipname=paste0(sub("Drosera_","",genome), "_", label)) %>%
       group_by(locus, tipname) %>% slice(1) %>% ungroup()
  all <- bind_rows(d, dio %>% select(locus,tip,tipname),
                   nep %>% select(locus,tip,tipname))
  loci <- all %>% count(locus) %>% filter(n >= MINGENE) %>% pull(locus)
  bind_rows(lapply(loci, function(l) {
    f <- file.path("DR/codon", paste0(l, ".fna"))
    if (!file.exists(f)) return(NULL)
    s <- readfa(f); a <- all[all$locus==l, ]
    a <- a[a$tip %in% names(s), ]
    if (nrow(a) < MINGENE) return(NULL)
    L <- nchar(s[a$tip[1]])
    if (L < 150 || L %% 3 != 0) return(NULL)
    tibble(locus=l, tipname=a$tipname, seq=unname(s[a$tip]), len=L)
  }))
}
cat("\n  building FULL set...\n"); flush.console()
FULL <- build(FALSE)
cat(sprintf("    loci: %d | total nt: %s\n", n_distinct(FULL$locus),
            format(sum(FULL$len[!duplicated(FULL$locus)]), big.mark=",")))
cat("  building STABLE-only set...\n"); flush.console()
STBL <- build(TRUE)
cat(sprintf("    loci: %d | total nt: %s\n", n_distinct(STBL$locus),
            format(sum(STBL$len[!duplicated(STBL$locus)]), big.mark=",")))

hr("3. WRITE SUPERMATRICES")
writephy <- function(D, subg, path) {
  keep <- c(paste0(sub("Drosera_","",DROS), "_", subg),
            paste0("Dionaea_", subg), "Nepenthes")
  d <- D %>% filter(tipname %in% keep)
  loci <- d %>% count(locus) %>% filter(n >= 4) %>% pull(locus)
  d <- d %>% filter(locus %in% loci)
  L <- d %>% distinct(locus, len) %>% arrange(locus)
  mat <- setNames(vapply(keep, function(t) {
    paste0(vapply(L$locus, function(l) {
      v <- d$seq[d$locus==l & d$tipname==t]
      if (length(v)) v[1] else strrep("-", L$len[L$locus==l])
    }, character(1)), collapse="")
  }, character(1)), keep)
  writeLines(c(sprintf(" %d %d", length(mat), nchar(mat[1])),
               paste0(names(mat), "  ", mat)), path)
  occ <- vapply(mat, function(x) 1 - nchar(gsub("[^-]","",x))/nchar(x), numeric(1))
  tibble(matrix=basename(path), tips=length(mat), sites=nchar(mat[1]),
         loci=nrow(L), min_occupancy=round(min(occ),3))
}
writeall <- function(D, path) {
  d <- D %>% filter(tipname %in% TIPS)
  loci <- d %>% count(locus) %>% filter(n >= 8) %>% pull(locus)
  d <- d %>% filter(locus %in% loci)
  L <- d %>% distinct(locus, len) %>% arrange(locus)
  mat <- setNames(vapply(TIPS, function(t) {
    paste0(vapply(L$locus, function(l) {
      v <- d$seq[d$locus==l & d$tipname==t]
      if (length(v)) v[1] else strrep("-", L$len[L$locus==l])
    }, character(1)), collapse="")
  }, character(1)), TIPS)
  writeLines(c(sprintf(" %d %d", length(mat), nchar(mat[1])),
               paste0(names(mat), "  ", mat)), path)
  occ <- vapply(mat, function(x) 1 - nchar(gsub("[^-]","",x))/nchar(x), numeric(1))
  cat("\n  per-tip occupancy in the 13-tip matrix:\n")
  print(round(sort(occ), 3))
  tibble(matrix=basename(path), tips=length(mat), sites=nchar(mat[1]),
         loci=nrow(L), min_occupancy=round(min(occ),3))
}
cat("\n  THE 13-TIP MATRIX is the one that answers DR03b's open question:\n")
cat("  it carries BOTH Dionaea subgenomes and BOTH subgenomes of every\n")
cat("  Drosera, so regia_A can be placed relative to Dionaea_A directly.\n")
cat("  The 7-tip A-only and B-only matrices are for the CONGRUENCE test:\n")
cat("  same species, two independent partitions, same topology expected\n")
cat("  under a single hybridisation.\n")
ST <- bind_rows(
  writeall(FULL, "DR/tree/all13_full.phy"),
  writeall(STBL, "DR/tree/all13_stable.phy"),
  writephy(FULL, "A", "DR/tree/A_full.phy"),
  writephy(FULL, "B", "DR/tree/B_full.phy"),
  writephy(STBL, "A", "DR/tree/A_stable.phy"),
  writephy(STBL, "B", "DR/tree/B_stable.phy"))
print(as.data.frame(ST), row.names=FALSE)
cat("\n  min_occupancy is the LEAST complete tip. Below ~0.3 that tip is\n")
cat("  mostly gaps and its placement will be poorly supported.\n")
write_csv(ST, "DR/out/DR05_matrix_stats.csv")

hr("4. PER-LOCUS ALIGNMENTS for the ASTRAL gene trees")
cat("  ASTRAL is consistent under the multi-species coalescent, so it handles\n")
cat("  the ILS that concatenation cannot. Needs one tree per locus.\n\n")
dir.create("DR/tree/genes", showWarnings=FALSE, recursive=TRUE)
unlink(list.files("DR/tree/genes", full.names=TRUE))
n <- 0
for (l in unique(FULL$locus)) {
  d <- FULL[FULL$locus==l, ]
  if (n_distinct(d$tipname) < 6) next
  writeLines(paste0(">", d$tipname, "\n", d$seq),
             file.path("DR/tree/genes", paste0(l, ".fna")))
  n <- n + 1
  if (n %% 200 == 0) { cat(sprintf("    %d gene alignments written\n", n)); flush.console() }
}
cat(sprintf("  %d per-locus alignments in DR/tree/genes/\n", n))
cat("\n  NOTE: these carry BOTH subgenomes, so each gene tree contains up to\n")
cat("  13 tips. ASTRAL is run on the A tips and the B tips separately.\n")

hr("DONE — next: bash DR/DR05_trees.sh")
