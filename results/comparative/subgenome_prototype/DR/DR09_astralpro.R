#!/usr/bin/env Rscript
# ============================================================================
# DR09 — SPECIES TREE FROM MULTI-COPY GENE TREES (ASTRAL-Pro)
#
# WHY
#   DR05 collapsed each species' A copies into ONE chimeric tip: where a
#   species had two A copies at a locus, one was picked arbitrarily. The
#   standing objection to regia_A's basal placement in the concatenated tree
#   has been that the chimera might be driving it. DR05c tested that by
#   restricting to loci where regia has exactly one A copy and got the same
#   topology -- but that only removed chimerism for REGIA, not for the others.
#
#   ASTRAL-Pro never collapses anything. It takes multi-copy gene trees and
#   models duplication and loss directly. If it ALSO places regia with
#   Dionaea, the chimerism objection is closed.
#
# WHAT CHANGES vs THE v5 INPUT
#   ASTRAL v5 needs one tip per SPECIES. ASTRAL-Pro needs one tip per GENE,
#   with duplicate species labels allowed, plus a gene->species map. So the
#   pruning we did for DR05 must NOT be applied here.
#
# THE SPECIES UNITS
#   Each subgenome is treated as its own "species": regia_A, regia_B, and so
#   on. That is the same convention as DR05 and it is what makes the tree
#   comparable. Within a unit, multiple gene copies are allowed and are what
#   ASTRAL-Pro is for.
#
# DATA
#   1,258 loci have at least one species with 2+ A copies; the median locus
#   carries 6 A copies pooled across species. ASTRAL-Pro uses all of it --
#   it does not require every species to be multi-copy at every locus, so the
#   usable set is ~1,258 rather than the 75 loci where all five are.
#
# IN   DR/tree/genetrees/*.treefile, DR/out/DR02_gene_delta.csv,
#      DR/out/DR02_segments.csv, DR/locus_meta.tsv,
#      fractionation_by_chrpair.csv
# OUT  DR/tree/pro/genetrees_pro.tre, DR/tree/pro/gene2species.map
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

hr("1. MAP EVERY GENE COPY TO A SUBGENOME UNIT")
LOC <- read_tsv("DR/locus_meta.tsv", show_col_types=FALSE) %>%
       group_by(tip) %>% slice(1) %>% ungroup()
P <- read_csv("fractionation_by_chrpair.csv", show_col_types=FALSE)
s1 <- P$retained_more == P$chrA
SIDE <- c(setNames(ifelse(s1,"A","B"), P$chrA), setNames(ifelse(s1,"B","A"), P$chrB))
dio <- LOC %>% filter(genome=="Dionaea_muscipula") %>%
  mutate(unit=paste0("Dionaea_", unname(SIDE[chr]))) %>%
  filter(!is.na(unname(SIDE[chr]))) %>% select(tip, unit)
nep <- LOC %>% filter(genome=="Nepenthes_gracilis") %>%
       transmute(tip, unit="Nepenthes")
G  <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE)
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE)
lab <- SS %>% filter(label %in% c("A","B")) %>%
       select(genome, region, chr, mb_lo, mb_hi, label)
dros <- G %>% mutate(mb=mid/1e6) %>%
  inner_join(lab, by=c("genome","region","chr"), relationship="many-to-many") %>%
  filter(mb>=mb_lo, mb<=mb_hi) %>% group_by(tip) %>% slice(1) %>% ungroup() %>%
  transmute(tip, unit=paste0(sub("Drosera_","",genome), "_", label))
MAP <- bind_rows(dio, nep, dros) %>% distinct(tip, .keep_all=TRUE)
cat(sprintf("  gene copies mapped to a unit: %d\n", nrow(MAP)))
print(as.data.frame(MAP %>% count(unit, name="copies") %>% arrange(unit)),
      row.names=FALSE)

hr("2. REBUILD GENE TREES WITH EVERY COPY KEPT")
cat("  DR05's gene trees were pruned to one tip per subgenome. ASTRAL-Pro\n")
cat("  wants them ALL -- duplicate labels are the point.\n\n")
fs <- list.files("DR/tree/pro/genetrees", pattern="[.]treefile$", full.names=TRUE)
cat(sprintf("  gene tree files: %d\n", length(fs)))
MAP$tip_iq <- gsub("@", "_", MAP$tip)
stopifnot(!any(duplicated(MAP$tip_iq)))
lk <- setNames(MAP$unit, MAP$tip_iq)
cat("  NOTE: IQ-TREE rewrites '@' as '_' in taxon names, so the map is keyed\n")
cat("  on the sanitised form. Example:", MAP$tip[1], "->", MAP$tip_iq[1], "\n\n")
keep <- 0L; multi <- 0L
con <- file("DR/tree/pro/genetrees_pro.tre", "w")
used <- character(0)
for (f in fs) {
  nw <- tryCatch(readLines(f, warn=FALSE)[1], error=function(e) NA)
  if (is.na(nw)) next
  tips <- regmatches(nw, gregexpr("[^(),:;]+(?=:)", nw, perl=TRUE))[[1]]
  tips <- trimws(tips)
  tips <- tips[tips %in% names(lk)]
  if (length(tips) < 4) next
  if (any(duplicated(lk[tips]))) multi <- multi + 1L
  writeLines(nw, con); keep <- keep + 1L
  used <- union(used, tips)
}
close(con)
cat(sprintf("  trees written: %d | of which multi-copy for some unit: %d\n",
            keep, multi))
write_tsv(MAP %>% filter(tip_iq %in% used) %>% select(tip_iq, unit),
          "DR/tree/pro/gene2species.map", col_names=FALSE)
cat(sprintf("  map written: %d gene->unit rows\n", sum(MAP$tip_iq %in% used)))
cat("\n  NOTE: DR05's gene trees kept ONE tip per subgenome, so multi-copy\n")
cat("  count above may be 0. If it is, section 3 explains the fix.\n")

hr("3. IS THE EXISTING GENE-TREE SET USABLE?")
if (multi == 0) {
  cat("  *** The existing gene trees were built from DR05's PRUNED alignments,\n")
  cat("  which already collapsed each subgenome to one copy. They cannot show\n")
  cat("  ASTRAL-Pro anything it does not already see. ***\n\n")
  cat("  The alignments in DR/codon/ still hold EVERY copy, so gene trees must\n")
  cat("  be rebuilt from those. Run:\n\n")
  cat("    ls DR/codon/*.fna | sed 's/.*\\///;s/[.]fna//' > /tmp/pro_loci.txt\n")
  cat("    mkdir -p DR/tree/pro/genetrees\n")
  cat("    parallel -j 64 --bar '/opt/share/software/bin/iqtree2 -s DR/codon/{}.fna \\\n")
  cat("      -m GTR+G -T 1 --prefix DR/tree/pro/genetrees/{} -redo >/dev/null 2>&1' \\\n")
  cat("      :::: /tmp/pro_loci.txt\n\n")
  cat("  then re-run this script, pointing section 2 at DR/tree/pro/genetrees\n")
} else {
  cat(sprintf("  %d trees carry multiple copies for at least one unit. Usable.\n", multi))
}
