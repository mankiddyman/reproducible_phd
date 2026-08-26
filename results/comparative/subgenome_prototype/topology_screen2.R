#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Cross-pair coherence screen using real OrthoFinder bitscores (Blast*.txt.gz, ofID keyed).
# CONFOUND: bitscore tracks rate as well as topology. Nepenthes sits outside both
# progenitors, so its skew is the pure rate/composition baseline. Topology signal
# = Drosera skew MINUS Nepenthes skew, measured per HOG where both are present.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
res <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE)
part <- res %>% filter(anchor == "Nepenthes_gracilis") %>%
  transmute(pair_lab, X = winner, Y = ifelse(winner == chrA, chrB, chrA))

sid <- readLines(file.path(GSD, "SpeciesIDs.txt"))
sid <- sid[nzchar(sid)]
spmap <- data.frame(idx = sub(":.*", "", sid),
                    sp  = sub("\\.(fa|faa|fasta)$", "", trimws(sub("^[0-9]+:", "", sid))),
                    stringsAsFactors = FALSE)
cat("=== SpeciesIDs ===\n"); print(spmap)
DIO <- spmap$idx[spmap$sp == "Dionaea_muscipula"]
stopifnot(length(DIO) == 1)

b <- bed %>% mutate(isRep = as.logical(isArrayRep)) %>% filter(is.na(isRep) | isRep)
cat(sprintf("ofID example: %s (gene %s)\n",
            b$ofID[b$genome == "Dionaea_muscipula"][1],
            b$id[b$genome == "Dionaea_muscipula"][1]))

xy <- bind_rows(transmute(part, chr = X, side = "X", pair_lab),
                transmute(part, chr = Y, side = "Y", pair_lab))
dio <- b %>% filter(genome == "Dionaea_muscipula") %>%
  inner_join(xy, by = "chr") %>% select(globHOG, ofID, side, pair_lab)

clean <- dio %>% group_by(globHOG) %>%
  filter(n() == 2, n_distinct(side) == 2, n_distinct(pair_lab) == 1) %>% ungroup() %>%
  pivot_wider(id_cols = c(globHOG, pair_lab), names_from = side, values_from = ofID)
cat(sprintf("clean 1X:1Y HOGs: %d\n", nrow(clean)))

tgt <- bind_rows(transmute(clean, ofID = X, globHOG, pair_lab, side = "X"),
                 transmute(clean, ofID = Y, globHOG, pair_lab, side = "Y"))
tf <- "topology/targets.txt"
writeLines(unique(tgt$ofID), tf)

## outgroup gene -> its HOG, so we can pair Drosera and Nepenthes within a HOG
og <- b %>% filter(genome != "Dionaea_muscipula") %>% select(ofID, globHOG, genome)

grab <- function(sp) {
  i <- spmap$idx[spmap$sp == sp]
  bf <- file.path(GSD, sprintf("Blast%s_%s.txt.gz", i, DIO))
  if (!file.exists(bf)) { message("missing ", bf); return(NULL) }
  cmd <- sprintf("zcat '%s' | awk -F'\\t' 'NR==FNR{t[$1]=1;next} ($2 in t){print $1\"\\t\"$2\"\\t\"$NF}' '%s' -",
                 bf, tf)
  x <- tryCatch(read.table(pipe(cmd), sep = "\t",
                           col.names = c("q", "s", "bits"), stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(x) || !nrow(x)) { message("no hits for ", sp); return(NULL) }
  x <- x %>% inner_join(select(og, ofID, globHOG), by = c("q" = "ofID")) %>%
    inner_join(select(tgt, ofID, side, pair_lab, globHOG), by = c("s" = "ofID", "globHOG"))
  bx <- x %>% filter(side == "X") %>% group_by(q, globHOG, pair_lab) %>%
    summarise(bX = max(bits), .groups = "drop")
  by <- x %>% filter(side == "Y") %>% group_by(q, globHOG, pair_lab) %>%
    summarise(bY = max(bits), .groups = "drop")
  inner_join(bx, by, by = c("q", "globHOG", "pair_lab")) %>%
    mutate(rel = (bX - bY) / pmax(bX, bY), sp = sp)
}

sps <- c(grep("^Drosera", spmap$sp, value = TRUE), "Nepenthes_gracilis")
v <- bind_rows(lapply(sps, grab))
if (!nrow(v)) stop("still no votes — paste the diagnostic output")
write_csv(v, "topology_votes.csv")

cat("\n=== marginal skew per species (NOT yet rate-corrected) ===\n")
sm <- v %>% group_by(sp) %>%
  summarise(n = n(), frac_X = mean(rel > 0), med_rel = median(rel), .groups = "drop")
print(as.data.frame(sm), digits = 3)

cat("\n=== per pair, fraction voting X ===\n")
pp <- v %>% group_by(sp, pair_lab) %>% summarise(f = mean(rel > 0), n = n(), .groups = "drop")
print(as.data.frame(pp %>% select(sp, pair_lab, f) %>%
                      pivot_wider(names_from = pair_lab, values_from = f)), digits = 3)

## ---- rate-corrected: per-HOG Drosera minus Nepenthes ----
nep <- v %>% filter(sp == "Nepenthes_gracilis") %>%
  group_by(globHOG, pair_lab) %>% summarise(rel_nep = median(rel), .groups = "drop")
dr <- v %>% filter(grepl("^Drosera", sp)) %>%
  group_by(globHOG, pair_lab) %>% summarise(rel_dros = median(rel), .groups = "drop")
pairw <- inner_join(dr, nep, by = c("globHOG", "pair_lab")) %>%
  mutate(delta = rel_dros - rel_nep)
cat(sprintf("\nHOGs with both Drosera and Nepenthes: %d\n", nrow(pairw)))

cat("\n=== rate-corrected topology signal (Drosera - Nepenthes, per HOG) ===\n")
tp <- pairw %>% group_by(pair_lab) %>%
  summarise(n = n(), med_delta = median(delta), frac_pos = mean(delta > 0),
            p = suppressWarnings(wilcox.test(delta)$p.value), .groups = "drop") %>%
  mutate(p_adj = p.adjust(p, "BH"))
print(as.data.frame(tp), digits = 3)
k <- sum(tp$med_delta > 0)
cat(sprintf("\nsame direction in %d/8 pairs | sign test p=%.4f\n",
            max(k, 8 - k), binom.test(max(k, 8 - k), 8, 0.5)$p.value))
cat(sprintf("pooled: median delta %.4f | Wilcoxon p=%.3g\n",
            median(pairw$delta), suppressWarnings(wilcox.test(pairw$delta)$p.value)))
write_csv(pairw, "topology_rate_corrected.csv")
