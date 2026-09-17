#!/usr/bin/env Rscript
# ============================================================================
# DR00 — BUILD THE LOCUS SET. No analysis, no filtering for a specific test.
#
# WHY A REBUILD
#   The existing wgd7 set (861 loci) required BOTH a Nepenthes anchor AND
#   Dionaea == 2 copies. Neither is needed by every downstream test:
#     - the ABC four-point (DR01) needs Nepenthes, NOT Dionaea
#     - delta labelling (DR02) needs Dionaea homeologs, NOT Nepenthes
#   Carrying both filters everywhere costs ~50% of the ABC loci and ~60% of
#   the labellable loci for no reason.
#
# THE UNION
#   SET N  every Nepenthes-anchored locus from synteny_ortho_hits.csv, with
#          ALL copies of every species. No Dionaea filter.       ~5,189
#   SET O  every Dionaea homeolog pair from the `og` column (one gene on each
#          side of one chromosome pair), plus all Drosera genes in the same
#          og group, plus Nepenthes if the og contains one.       ~2,004
#          ~1,005 of these already appear in SET N.
#   UNION  ~6,188 loci, of which 867 are already built.
#
# WHY `og` AND NOT `globHOG`
#   GENESPACE's `og` is the SYNTENY-CONSTRAINED orthogroup (built by
#   build_synOGs); globOG/globHOG are genome-wide OrthoFinder output.
#   Lovell et al. 2022 eLife show that for two meso-tetraploid cottons sharing
#   a WGD that predated speciation, single-copy orthogroups were 35.6% under
#   default OrthoFinder but 85.7% under GENESPACE syntenic orthogroups. Same
#   situation here: og gives 2,004 clean Dionaea homeolog pairs vs globHOG's
#   1,245, and 96% of two-sided og groups have exactly one gene per side
#   against 63% for globHOG.
#
# NAMING
#   tip  = genome@gene_id, matching wgd/all_pep.fa and cds/all_cds_tagged.fa
#   locus id = the Nepenthes gene where one exists, else "og<number>"
#   Written to DR/, NOT wgd7/ — the old set stays intact for comparison.
#
# IN   synteny_ortho_hits.csv, ../genespace/results/combBed.txt,
#      fractionation_by_chrpair.csv, wgd/all_pep.fa
# OUT  DR/ids/<locus>.ids, DR/locus_meta.tsv, DR/anchorlist.txt,
#      DR/out/DR00_locus_summary.csv
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE)
DR <- "DR"; GSD <- file.path(dirname(BASE),"genespace","results")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")
ALL <- c("Nepenthes_gracilis","Dionaea_muscipula", DROS)

hr("1. SET N — Nepenthes-anchored, ALL copies, no Dionaea filter")
h <- read_csv("synteny_ortho_hits.csv", show_col_types=FALSE) %>%
     filter(sp %in% c("Dionaea_muscipula", DROS))
cat(sprintf("  anchors: %d | rows: %d\n", n_distinct(h$nep_gene), nrow(h)))
cat("\n  copies per anchor by species (the Dio==2 filter is NOT applied):\n")
print(as.data.frame(h %>% count(nep_gene, sp, name="k") %>% group_by(sp) %>%
  summarise(anchors=n(), median_k=median(k), mean_k=round(mean(k),2),
            p95=quantile(k,.95), max_k=max(k), .groups="drop")))
setN <- bind_rows(
  h %>% distinct(nep_gene, nep_chr) %>%
    transmute(locus=nep_gene, tip=paste0("Nepenthes_gracilis@", nep_gene),
              gene=nep_gene, genome="Nepenthes_gracilis", chr=nep_chr),
  h %>% transmute(locus=nep_gene, tip=paste0(sp,"@",sp_gene),
                  gene=sp_gene, genome=sp, chr=sp_chr)) %>%
  mutate(set="N", has_nep=TRUE)
cat(sprintf("\n  SET N: %d loci, %d tips\n", n_distinct(setN$locus), nrow(setN)))

hr("2. SET O — Dionaea homeolog pairs from the syntenic orthogroup column")
bed <- read.table(file.path(GSD,"combBed.txt"), header=TRUE, sep="\t",
                  quote="", comment.char="", stringsAsFactors=FALSE) %>%
       mutate(isRep=as.logical(isArrayRep)) %>% filter(is.na(isRep)|isRep)
PAIRS <- read_csv("fractionation_by_chrpair.csv", show_col_types=FALSE)
SIDE <- c(setNames(rep("A",8),PAIRS$chrA), setNames(rep("B",8),PAIRS$chrB))
PAIR <- c(setNames(PAIRS$exp_pair,PAIRS$chrA), setNames(PAIRS$exp_pair,PAIRS$chrB))
hp <- bed %>% filter(genome=="Dionaea_muscipula") %>%
  mutate(side=unname(SIDE[chr]), pr=unname(PAIR[chr])) %>% filter(!is.na(side)) %>%
  group_by(og) %>%
  summarise(nA=sum(side=="A"), nB=sum(side=="B"), np=n_distinct(pr), .groups="drop") %>%
  filter(nA==1, nB==1, np==1)
cat(sprintf("  Dionaea homeolog pairs (og): %d\n", nrow(hp)))
setO_all <- bed %>% filter(og %in% hp$og, genome %in% ALL) %>%
  transmute(og, tip=paste0(genome,"@",id), gene=id, genome, chr)
nep_in_og <- setO_all %>% filter(genome=="Nepenthes_gracilis") %>%
             distinct(og, nep_gene=gene)
cat(sprintf("  og pairs containing a Nepenthes gene: %d (%.0f%%)\n",
            nrow(nep_in_og), 100*nrow(nep_in_og)/nrow(hp)))
cat("  (those without one are DELTA-ONLY -- no outgroup, so no four-point.)\n")
setO <- setO_all %>% left_join(nep_in_og, by="og") %>%
  mutate(locus=ifelse(!is.na(nep_gene), nep_gene, paste0("og", og)),
         set="O", has_nep=!is.na(nep_gene)) %>% select(-og, -nep_gene)

hr("3. UNION")
LOC <- bind_rows(setN, setO) %>% distinct(locus, tip, .keep_all=TRUE) %>%
  group_by(locus) %>%
  mutate(has_nep=any(genome=="Nepenthes_gracilis"),
         set=ifelse(n_distinct(set)>1, "both", set[1])) %>% ungroup()
cat(sprintf("  loci: %d | tips: %d\n", n_distinct(LOC$locus), nrow(LOC)))
print(as.data.frame(LOC %>% distinct(locus, set, has_nep) %>%
  count(set, has_nep, name="loci")))
cat("\n  tips per locus:\n"); print(summary(as.integer(table(LOC$locus))))
tl <- as.integer(table(LOC$locus))
cat(sprintf("  loci with >30 tips: %d | >40: %d  (mafft cost scales ~n^2)\n",
            sum(tl>30), sum(tl>40)))
cat("\n  loci usable per test:\n")
cap <- LOC %>% group_by(locus) %>%
  summarise(nep=any(genome=="Nepenthes_gracilis"),
            dio=sum(genome=="Dionaea_muscipula"), .groups="drop")
cat(sprintf("    four-point (needs Nepenthes)          : %d\n", sum(cap$nep)))
cat(sprintf("    delta (needs exactly 2 Dionaea copies): %d\n", sum(cap$dio==2)))
cat(sprintf("    both                                  : %d\n", sum(cap$nep & cap$dio==2)))

hr("4. VALIDATE — do all tip names exist in the peptide fasta?")
hdr <- sub("^>", "", grep("^>", readLines("wgd/all_pep.fa"), value=TRUE))
hdr <- sub(" .*$", "", hdr)
miss <- setdiff(unique(LOC$tip), hdr)
cat(sprintf("  tips: %d | present in all_pep.fa: %d | MISSING: %d\n",
            n_distinct(LOC$tip), n_distinct(LOC$tip)-length(miss), length(miss)))
if (length(miss)) { cat("  first missing:\n"); print(head(miss, 10))
  cat("  *** STOPPING — name construction is wrong. ***\n"); quit(status=1) }
cat("  all tip names resolve.\n")

hr("5. WRITE")
already <- sub("\\.fna$","", list.files("wgd7/codon_gap", pattern="\\.fna$"))
todo <- setdiff(unique(LOC$locus), already)
cat(sprintf("  loci already built in wgd7: %d\n", length(intersect(unique(LOC$locus), already))))
cat(sprintf("  loci to build: %d\n", length(todo)))
cat(sprintf("  estimated mafft at 0.50 s/locus, -j 40: %.1f min\n", length(todo)*0.50/40/60))
dir.create(file.path(DR,"ids"), recursive=TRUE, showWarnings=FALSE)
for (a in unique(LOC$locus))
  writeLines(LOC$tip[LOC$locus==a], file.path(DR,"ids", paste0(a,".ids")))
write_tsv(LOC, file.path(DR,"locus_meta.tsv"))
writeLines(sort(unique(LOC$locus)), file.path(DR,"anchorlist.txt"))
writeLines(sort(todo), file.path(DR,"buildlist.txt"))
write_csv(LOC %>% distinct(locus, set, has_nep) %>%
  left_join(cap, by="locus") %>% rename(n_dionaea=dio),
  file.path(DR,"out","DR00_locus_summary.csv"))
cat(sprintf("\n  wrote %d .ids files, anchorlist.txt, buildlist.txt\n", n_distinct(LOC$locus)))
