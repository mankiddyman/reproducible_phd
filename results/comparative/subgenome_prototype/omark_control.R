#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# OMArk classifies genes against a reference orthology DB, independent of
# GENESPACE/OrthoFinder/synteny. Two questions:
#   1. Is Y enriched for fragmented/erroneous gene models? (annotation quality)
#   2. Did the OMArk cleaning step remove genes asymmetrically? (manufactured signal)
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

AD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/Dionaea_muscipula/annotation"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

res  <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE)
part <- res %>% filter(anchor == "Nepenthes_gracilis") %>%
  transmute(pair_lab, X = winner, Y = ifelse(winner == chrA, chrB, chrA))
xy <- bind_rows(transmute(part, chr = X, side = "X", pair_lab),
                transmute(part, chr = Y, side = "Y", pair_lab))

gff_genes <- function(p) {
  x <- read.table(p, sep = "\t", quote = "", comment.char = "#",
                  stringsAsFactors = FALSE, fill = TRUE)
  x %>% filter(V3 == "gene") %>%
    transmute(chr = V1, gene = sub(";.*", "", sub(".*ID=", "", V9)))
}
raw   <- gff_genes(file.path(AD, "annevo/Dionaea_muscipula.annevo.gff3"))
clean <- gff_genes(file.path(AD, "annevo/Dionaea_muscipula.annevo.omark_clean.gff3"))
cat(sprintf("genes: raw %d | omark_clean %d | removed %d\n",
            nrow(raw), nrow(clean), nrow(raw) - nrow(clean)))

## ---- did OMArk cleaning hit X and Y differently? ----
rem <- raw %>% mutate(kept = gene %in% clean$gene) %>% inner_join(xy, by = "chr") %>%
  group_by(pair_lab, side) %>%
  summarise(n_raw = n(), n_removed = sum(!kept), frac_removed = mean(!kept), .groups = "drop") %>%
  pivot_wider(names_from = side, values_from = c(n_raw, n_removed, frac_removed))
cat("\n=== OMArk removals, X vs Y ===\n"); print(as.data.frame(rem), digits = 3)
if (all(c("frac_removed_X","frac_removed_Y") %in% names(rem)))
  cat(sprintf("Y removed more in %d/8 | paired Wilcoxon p=%.3f  (Y higher => cleaning biased against Y)\n",
              sum(rem$frac_removed_Y > rem$frac_removed_X, na.rm = TRUE),
              suppressWarnings(wilcox.test(rem$frac_removed_X, rem$frac_removed_Y,
                                           paired = TRUE)$p.value)))

## ---- OMArk class composition per side ----
gcf <- file.path(AD, "annevo/Dionaea_muscipula.annevo.gene_class.tsv")
gc  <- read.table(gcf, sep = "\t", header = TRUE, quote = "", comment.char = "",
                  stringsAsFactors = FALSE)
names(gc)[1] <- "gene"
cls <- if ("reason" %in% names(gc)) "reason" else names(gc)[ncol(gc)]
cat(sprintf("\nusing OMArk class column: %s\n", cls))

j <- raw %>% inner_join(gc, by = "gene") %>% inner_join(xy, by = "chr")
tabl <- j %>% count(side, .data[[cls]]) %>%
  group_by(side) %>% mutate(frac = n/sum(n)) %>% ungroup()
cat("\n=== OMArk class composition, X vs Y ===\n")
print(as.data.frame(tabl %>% select(side, class = all_of(cls), n, frac)), digits = 3)

per <- j %>% count(pair_lab, side, .data[[cls]]) %>%
  group_by(pair_lab, side) %>% mutate(frac = n/sum(n)) %>% ungroup()
for (k in unique(per[[cls]])) {
  w <- per %>% filter(.data[[cls]] == k) %>% select(pair_lab, side, frac) %>%
    pivot_wider(names_from = side, values_from = frac)
  if (!all(c("X","Y") %in% names(w)) || sum(complete.cases(w)) < 5) next
  cat(sprintf("%-22s Y>X in %d/8 | paired p=%.3f\n", k,
              sum(w$Y > w$X, na.rm = TRUE),
              suppressWarnings(wilcox.test(w$X, w$Y, paired = TRUE)$p.value)))
}
write_csv(rem, "omark_removals.csv"); write_csv(per, "omark_classes.csv")
