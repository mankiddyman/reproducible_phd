#!/usr/bin/env Rscript
# ============================================================================
# DR12 — HOW MUCH OF THE B SUBGENOME IS INTROGRESSED?
#
# DR11 found D = 0.296 (p = 6.8e-14) in B and nothing in A, by counting gene
# tree topologies. Two limitations of that:
#   1. its p-value treats gene trees as independent
#   2. D measures ASYMMETRY, not admixture proportion. A small amount of gene
#      flow can produce a large D.
#
# THIS SCRIPT
#   (a) recomputes D from SITE PATTERNS with a block jackknife over loci,
#       which is the standard calibration and accounts for linkage
#   (b) reports what D can and cannot be converted into
#
# ONE THING DR11 ALREADY TELLS US
#   D per core species: 0.280, 0.297, 0.309, 0.297 -- near identical. If the
#   gene flow came from one sampled lineage, that species would stand out.
#   It does not, so the source is the COMMON ANCESTOR of core Drosera or a
#   lineage close to it. That rules out the standard f4-ratio, which needs an
#   unadmixed sister to the source.
#
# THE CONFIGURATION
#   (((P1, P2), P3), O)  with  P1 = Dionaea_B, P2 = regia_B,
#                              P3 = a core Drosera B, O = Nepenthes
#   ABBA = P2 and P3 share a derived state   -> regia allied with core
#   BABA = P1 and P3 share a derived state   -> Dionaea allied with core
#   D = (ABBA - BABA)/(ABBA + BABA);  0 under pure ILS
#
# IN   DR/tree/B_full.phy, DR/tree/A_full.phy (control)
# OUT  DR/out/DR12_dstat.csv
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
CORE <- c("binata","paradoxa","scorpioides","capensis")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

readphy <- function(f) {
  x <- readLines(f, warn=FALSE); x <- x[nzchar(trimws(x))][-1]
  nm <- sub("\\s.*$", "", trimws(x))
  sq <- toupper(gsub("^\\S+\\s+", "", trimws(x)))
  setNames(strsplit(sq, ""), nm)
}
dstat <- function(S, p1, p2, p3, o, nblock=50) {
  if (!all(c(p1,p2,p3,o) %in% names(S))) return(NULL)
  a <- S[[p1]]; b <- S[[p2]]; c3 <- S[[p3]]; oo <- S[[o]]
  ok <- a %in% c("A","C","G","T") & b %in% c("A","C","G","T") &
        c3 %in% c("A","C","G","T") & oo %in% c("A","C","G","T")
  a <- a[ok]; b <- b[ok]; c3 <- c3[ok]; oo <- oo[ok]
  abba <- (a == oo) & (b == c3) & (b != oo)
  baba <- (b == oo) & (a == c3) & (a != oo)
  n <- length(a); blk <- cut(seq_len(n), nblock, labels=FALSE)
  A <- sum(abba); B <- sum(baba)
  if (A + B < 50) return(NULL)
  Dall <- (A - B)/(A + B)
  jk <- vapply(seq_len(nblock), function(k) {
    ai <- sum(abba[blk != k]); bi <- sum(baba[blk != k])
    if (ai + bi == 0) NA_real_ else (ai - bi)/(ai + bi)
  }, numeric(1))
  jk <- jk[!is.na(jk)]; m <- length(jk)
  se <- sqrt((m-1)/m * sum((jk - mean(jk))^2))
  tibble(P1=p1, P2=p2, P3=p3, sites=n, ABBA=A, BABA=B,
         D=Dall, SE=se, Z=Dall/se,
         p=2*pnorm(-abs(Dall/se)))
}

hr("1. SITE-PATTERN D, B SUBGENOME")
cat("  P1 = Dionaea_B, P2 = regia_B, P3 = each core species, O = Nepenthes\n")
cat("  D > 0 means regia_B shares more derived alleles with core Drosera\n")
cat("  than Dionaea_B does. Block jackknife over 50 blocks.\n\n")
SB <- readphy("DR/tree/B_full.phy")
cat("  taxa:", paste(names(SB), collapse=", "), "\n")
cat("  sites:", length(SB[[1]]), "\n\n")
RB <- bind_rows(lapply(CORE, function(cs)
  dstat(SB, "Dionaea_B", "regia_B", paste0(cs,"_B"), "Nepenthes"))) %>%
  mutate(subgenome="B")
print(as.data.frame(RB %>% transmute(P3, sites, ABBA, BABA,
  D=round(D,4), SE=round(SE,4), Z=round(Z,2), p=signif(p,3))), row.names=FALSE)

hr("2. THE CONTROL — A SUBGENOME")
cat("  Same test, same code. Should give D near 0.\n\n")
SA <- readphy("DR/tree/A_full.phy")
RA <- bind_rows(lapply(CORE, function(cs)
  dstat(SA, "Dionaea_A", "regia_A", paste0(cs,"_A"), "Nepenthes"))) %>%
  mutate(subgenome="A")
print(as.data.frame(RA %>% transmute(P3, sites, ABBA, BABA,
  D=round(D,4), SE=round(SE,4), Z=round(Z,2), p=signif(p,3))), row.names=FALSE)
write_csv(bind_rows(RB, RA), "DR/out/DR12_dstat.csv")

hr("3. DOES IT AGREE WITH THE TOPOLOGY-BASED D?")
cat(sprintf("  site-pattern D, B: %.3f (mean over four core species)\n",
            mean(RB$D)))
cat("  topology D,      B: 0.296\n")
cat(sprintf("  site-pattern D, A: %.3f\n", mean(RA$D)))
cat("  topology D,      A: -0.014\n")
cat("\n  Agreement means both are measuring the same asymmetry, and the\n")
cat("  jackknife Z is the properly calibrated significance.\n")

hr("4. WHAT D CANNOT TELL YOU")
cat("  D is the ASYMMETRY, not the admixture proportion. A small amount of\n")
cat("  gene flow between closely related lineages produces a large D.\n\n")
cat("  For a proportion you need an f4-RATIO, which requires an unadmixed\n")
cat("  sister to the SOURCE. DR11 showed all four core species give nearly\n")
cat("  identical D (0.280-0.309), so the source is their common ancestor and\n")
cat("  no sampled species can serve that role. The f4-ratio is unavailable.\n\n")
cat("  HyDe estimates gamma (inheritance probability) from phylogenetic\n")
cat("  invariants and does NOT need that sister. Run it next:\n\n")
cat("    python3 -m phyde.scripts.run_hyde -i DR/tree/B_full.phy \\\n")
cat("      -m DR/out/hyde_map.txt -o Nepenthes -n 7 -t 7 -s <nsites> \\\n")
cat("      --prefix DR/out/hyde_B\n\n")
m <- tibble(ind=names(SB), taxon=names(SB))
write_tsv(m, "DR/out/hyde_map.txt", col_names=FALSE)
cat(sprintf("  wrote DR/out/hyde_map.txt (%d taxa, one individual each)\n", nrow(m)))
cat(sprintf("  nsites for B = %d, for A = %d\n", length(SB[[1]]), length(SA[[1]])))
