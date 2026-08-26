#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Gene content from the RAW annotation. No orthogroups, no synteny.
# Separates real gene-count asymmetry from amplification by HOG assignment.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

AD  <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/Dionaea_muscipula/annotation"
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

res  <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE)
part <- res %>% filter(anchor == "Nepenthes_gracilis") %>%
  transmute(pair_lab, X = winner, Y = ifelse(winner == chrA, chrB, chrA))
xy <- bind_rows(transmute(part, chr = X, side = "X", pair_lab),
                transmute(part, chr = Y, side = "Y", pair_lab))

gff <- read.table(file.path(AD, "annevo/Dionaea_muscipula.annevo.gff3"), sep = "\t",
                  quote = "", comment.char = "#", fill = TRUE, stringsAsFactors = FALSE) %>%
  filter(V3 == "gene") %>% transmute(chr = V1, gene = sub(";.*", "", sub(".*ID=", "", V9)))

gc <- read.table(file.path(AD, "annevo/Dionaea_muscipula.annevo.gene_class.tsv"),
                 sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
names(gc)[1] <- "gene"

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
norm <- function(x) sub("\\.(t|p)?[0-9]+$", "", sub("[-_]?(mRNA|RNA)[-_]?[0-9]+$", "", x))
hogset <- bed %>% filter(genome == "Dionaea_muscipula") %>%
  mutate(isRep = as.logical(isArrayRep)) %>% filter(is.na(isRep) | isRep) %>% pull(id)
hogset <- norm(hogset)
cat(sprintf("normalised combBed id example: %s\n", hogset[1]))
cat(sprintf("id match: gff->gene_class %.3f | gff->combBed %.3f\n",
            mean(gff$gene %in% gc$gene), mean(hogset %in% gff$gene)))

g <- gff %>% left_join(gc, by = "gene") %>% inner_join(xy, by = "chr") %>%
  mutate(in_hog = gene %in% hogset,
         keep_omark = decision == "keep",
         consistent = reason == "Consistent")

sm <- g %>% group_by(pair_lab, side) %>%
  summarise(n_all = n(), n_keep = sum(keep_omark, na.rm = TRUE),
            n_cons = sum(consistent, na.rm = TRUE),
            n_hog = sum(in_hog), hog_rate = mean(in_hog), .groups = "drop")

w <- sm %>% pivot_wider(names_from = side,
                        values_from = c(n_all, n_keep, n_cons, n_hog, hog_rate)) %>%
  mutate(r_all = n_all_X/n_all_Y, r_keep = n_keep_X/n_keep_Y,
         r_cons = n_cons_X/n_cons_Y, r_hog = n_hog_X/n_hog_Y)

cat("\n=== gene counts per side, raw annotation ===\n")
print(as.data.frame(select(w, pair_lab, n_all_X, n_all_Y, n_cons_X, n_cons_Y,
                           n_hog_X, n_hog_Y)))
cat("\n=== X/Y ratio under progressively stricter filters ===\n")
print(as.data.frame(select(w, pair_lab, r_all, r_keep, r_cons, r_hog)), digits = 3)
for (v in c("r_all","r_keep","r_cons","r_hog")) {
  x <- w[[v]][is.finite(w[[v]])]
  if (!length(x)) { cat(sprintf("%-7s no finite values\n", v)); next }
  k <- sum(x > 1); n <- length(x)
  cat(sprintf("%-7s X>Y in %d/%d | mean %.3f | sign-test p=%.4f\n", v, k, n, mean(x),
              binom.test(max(k, n - k), n, 0.5)$p.value))
}

cat("\n=== orthogroup-assignment rate (the amplification step) ===\n")
print(as.data.frame(select(w, pair_lab, hog_rate_X, hog_rate_Y)), digits = 3)
cat(sprintf("X>Y in %d/8 | paired Wilcoxon p=%.3f | mean X %.3f vs Y %.3f\n",
            sum(w$hog_rate_X > w$hog_rate_Y),
            suppressWarnings(wilcox.test(w$hog_rate_X, w$hog_rate_Y, paired = TRUE)$p.value),
            mean(w$hog_rate_X), mean(w$hog_rate_Y)))
cat(sprintf("\namplification: raw %.3f -> HOG %.3f (%.0f%% of excess added by HOG step)\n",
            mean(w$r_all), mean(w$r_hog),
            100*(mean(w$r_hog) - mean(w$r_all))/(mean(w$r_hog) - 1)))
write_csv(w, "raw_gene_content.csv")
