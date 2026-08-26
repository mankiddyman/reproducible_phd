#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
DIO <- "Dionaea_muscipula"
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")

h <- read_csv("synteny_ortho_hits.csv", show_col_types=FALSE)
hr("0. what is in synteny_ortho_hits.csv")
cat("  columns:", paste(names(h), collapse=", "), "\n")
cat(sprintf("  rows: %d   anchors: %d   species: %s\n",
            nrow(h), n_distinct(h$nep_gene), paste(sort(unique(h$sp)), collapse=" ")))

h <- filter(h, sp %in% c(DIO, DROS))
dio2 <- h %>% filter(sp==DIO) %>% count(nep_gene, name="k") %>%
        filter(k==2) %>% pull(nep_gene)
cat(sprintf("\n  anchors total          : %d\n", n_distinct(h$nep_gene)))
cat(sprintf("  anchors with Dio == 2  : %d  (%.1f%%)\n",
            length(dio2), 100*length(dio2)/n_distinct(h$nep_gene)))
cat("  ^ this is the wgd7 filter, and it is the binding constraint\n")

hr("1. usable loci for the within-species four-point test")
cat("  Needs: >=3 copies on >=3 distinct chromosomes. Dionaea is NOT required.\n\n")
usable <- function(d) d %>% group_by(sp, nep_gene) %>%
  summarise(k=n(), nchr=n_distinct(sp_chr), .groups="drop") %>%
  filter(k>=3, nchr>=3)
now <- usable(h %>% filter(nep_gene %in% dio2)) %>% count(sp, name="current")
all <- usable(h) %>% count(sp, name="available")
cmp <- full_join(all, now, by="sp") %>%
  mutate(current=coalesce(current,0L), gain=available-current,
         fold=round(available/pmax(1,current),1)) %>% arrange(desc(gain))
print(as.data.frame(cmp), row.names=FALSE)
cat("\n  'current' should match section 0 of script 24:\n")
cat("  binata 24, scorpioides 16, paradoxa 31, regia 139, capensis 105\n")

hr("2. is the extra material biased?")
cat("  Loci with Dio != 2 are loci where Dionaea fractionated. If fractionation\n")
cat("  is correlated across species, the new loci may differ systematically.\n\n")
cop <- h %>% filter(sp != DIO) %>% group_by(sp, nep_gene) %>%
  summarise(k=n(), nchr=n_distinct(sp_chr), .groups="drop") %>%
  mutate(set=ifelse(nep_gene %in% dio2, "Dio==2", "Dio!=2"))
print(as.data.frame(cop %>% group_by(sp, set) %>%
  summarise(anchors=n(), mean_k=round(mean(k),2), frac_k1=round(mean(k==1),3),
            frac_3chr=round(mean(nchr>=3),3), .groups="drop") %>%
  arrange(sp, set)), row.names=FALSE)
cat("\n  Similar frac_k1 and frac_3chr across the two sets = the extra loci are\n")
cat("  comparable and the extension is safe. Large differences = flag it.\n")

hr("3. power estimate for the AAB vs ABC test")
cat("  Script 24 section 2 needs >=8 triplets per region to test uniformity.\n")
cat("  Region here = Nepenthes _dom chromosome.\n\n")
reg <- h %>% filter(sp != DIO) %>% group_by(sp, nep_gene, nep_chr) %>%
  summarise(k=n(), nchr=n_distinct(sp_chr), .groups="drop") %>%
  filter(k>=3, nchr>=3) %>% count(sp, nep_chr, name="loci")
print(as.data.frame(reg %>% group_by(sp) %>%
  summarise(regions_total=n(), regions_ge8=sum(loci>=8),
            regions_ge15=sum(loci>=15), max_loci=max(loci), .groups="drop") %>%
  arrange(desc(regions_ge8))), row.names=FALSE)
cat("\n  regia currently has 3 regions at p_adj < 0.05 from 6 testable regions.\n")
cat("  A species needs several regions at >=8 for the same test to have teeth.\n")

hr("4. build cost")
need <- usable(h) %>% count(sp, name="loci")
cat("  One alignment per species-locus: that species' copies + 1 Nepenthes tip.\n")
cat("  4-6 tips each, so far cheaper per alignment than wgd7.\n\n")
print(as.data.frame(need), row.names=FALSE)
cat(sprintf("\n  total alignments to build: %d\n", sum(need$loci)))
cat("  pipeline: seqkit grep -> mafft --localpair -> pal2nal (gapped) -> yn00\n")
cat("  inputs already on disk: wgd/all_pep.fa, cds/all_cds_tagged.fa\n")
for (f in c("wgd/all_pep.fa","cds/all_cds_tagged.fa","all_pep.faa","all_cds.fna"))
  cat(sprintf("    %-28s %s\n", f, ifelse(file.exists(f), "present", "MISSING")))
