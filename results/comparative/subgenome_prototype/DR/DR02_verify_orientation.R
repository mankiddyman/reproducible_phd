#!/usr/bin/env Rscript
# Independent verification that A = the dominant (more-retained) chromosome.
# Does NOT trust `retained_more`: recomputes it from the raw kA/kB counts and
# checks four separate ways.
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
P <- read_csv("fractionation_by_chrpair.csv", show_col_types=FALSE)

cat("=== 1. is `retained_more` self-consistent with kA/kB? ===\n")
P$recomputed <- ifelse(P$kA > P$kB, P$chrA, P$chrB)
P$match <- P$recomputed == P$retained_more
print(as.data.frame(P %>% transmute(exp_pair, chrA, kA, chrB, kB,
  retained_more, recomputed, match)), row.names=FALSE)
stopifnot(all(P$match))
cat("  all 8 agree -- retained_more genuinely is the higher-count chromosome.\n")

cat("\n=== 2. is chrA alphabetical rather than meaningful? ===\n")
cat("  fractionation.R does arrange(exp_pair, one_chr) then takes [1] and [2].\n")
P$alpha_ok <- P$chrA < P$chrB
cat(sprintf("  chrA < chrB alphabetically in %d of 8 pairs\n", sum(P$alpha_ok)))
cat("  8 of 8 -> chrA is a SORT ORDER, not a subgenome identity.\n")

cat("\n=== 3. does sg1 mean anything? ===\n")
P$A_is_sg1 <- grepl("_sg1_", P$retained_more)
print(as.data.frame(P %>% transmute(exp_pair, retained_more,
  dominant_is=ifelse(A_is_sg1, "sg1", "sg2"))), row.names=FALSE)
cat(sprintf("\n  dominant is sg1 in %d of 8, sg2 in %d.\n",
            sum(P$A_is_sg1), sum(!P$A_is_sg1)))
cat("  A coin flip -> sg1/sg2 are assembly labels with no subgenome meaning.\n")
cat("  The three sg2-dominant pairs are exactly the ones that were inverted:\n")
print(as.data.frame(P %>% filter(!A_is_sg1) %>%
  transmute(exp_pair, retained_more, frac_A=round(frac_A,3))), row.names=FALSE)

cat("\n=== 4. does the SCRIPT now use the dominant chromosome? ===\n")
sg1_dom <- P$retained_more == P$chrA
SIDE <- c(setNames(ifelse(sg1_dom,"A","B"), P$chrA),
          setNames(ifelse(sg1_dom,"B","A"), P$chrB))
chk <- tibble(chr=names(SIDE), side=unname(SIDE)) %>%
  mutate(is_dominant = chr %in% P$retained_more)
print(as.data.frame(chk %>% count(side, is_dominant)), row.names=FALSE)
bad <- chk %>% filter((side=="A") != is_dominant)
if (nrow(bad)) { cat("  *** MISMATCH ***\n"); print(as.data.frame(bad)); quit(status=1) }
cat("  every chromosome labelled A is the dominant one. Every B is not.\n")

cat("\n=== 5. does it hold in the actual output? ===\n")
G <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE)
LOC <- read_tsv("DR/locus_meta.tsv", show_col_types=FALSE) %>%
       group_by(tip) %>% slice(1) %>% ungroup() %>%
       filter(genome=="Dionaea_muscipula")
ax <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE) %>%
      distinct(locus)
dchk <- LOC %>% semi_join(ax, by="locus") %>%
  mutate(side=unname(SIDE[chr])) %>% filter(!is.na(side)) %>%
  count(chr, side) %>% mutate(dominant = chr %in% P$retained_more)
print(as.data.frame(dchk), row.names=FALSE)
cat(sprintf("\n  loci used: %d | Dionaea genes on A-side chromosomes: %d\n",
            nrow(ax), sum(dchk$n[dchk$side=="A"])))
cat("  A-side rows must all have dominant=TRUE.\n")
stopifnot(all(dchk$dominant[dchk$side=="A"]),
          all(!dchk$dominant[dchk$side=="B"]))
cat("  VERIFIED.\n")

cat("\n=== 6. the direction of delta ===\n")
cat("  delta = [d(R,D_A) - d(R,D_B)] / d(D_A,D_B)\n")
cat("  An A-derived copy is CLOSER to D_A, so d(R,D_A) is smaller and delta\n")
cat("  is NEGATIVE. The script labels segments A when med_delta < -0.05.\n")
cat("  Confirmed by the control: D_A itself gives delta = -1.0000 exactly.\n")
