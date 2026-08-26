#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Chain HSPs per (query, target window), score query coverage, compare
# "lost" genes against the positive control (homolog present + annotated)
# and the negative control (non-homeologous chromosome).
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

jb <- read_csv("blast/jobs.csv", show_col_types = FALSE)
CN <- c("qseqid","qlen","sseqid","pident","length","qstart","qend","sstart","send","evalue","bits")

cov_union <- function(s, e) { o <- order(s); s <- s[o]; e <- e[o]
  cs <- s[1]; ce <- e[1]; tot <- 0
  if (length(s) > 1) for (i in 2:length(s)) {
    if (s[i] <= ce + 1) ce <- max(ce, e[i]) else { tot <- tot + ce - cs + 1; cs <- s[i]; ce <- e[i] } }
  tot + ce - cs + 1 }

score <- function(tag) {
  f <- file.path("blast/out", paste0(tag, ".tsv"))
  if (!file.exists(f) || file.size(f) == 0) return(NULL)
  h <- read.table(f, col.names = CN, stringsAsFactors = FALSE)
  h %>% mutate(win = floor(pmin(sstart, send) / 1e5)) %>%
    group_by(qseqid, qlen, sseqid, win) %>%
    summarise(cov = cov_union(qstart, qend)/qlen[1],
              pid = weighted.mean(pident, length), .groups = "drop") %>%
    group_by(qseqid) %>% slice_max(cov, n = 1, with_ties = FALSE) %>% ungroup() %>%
    mutate(tag = tag)
}

hits <- bind_rows(lapply(jb$tag, score))
if (!nrow(hits)) stop("no BLAST output yet")

full <- jb %>% rowwise() %>%
  mutate(ids = list(readLines(file.path("blast/prot", paste0(tag, ".ids"))))) %>%
  unnest(ids) %>% rename(qseqid = ids) %>%
  left_join(select(hits, qseqid, tag, cov, pid), by = c("qseqid","tag")) %>%
  mutate(cov = replace_na(cov, 0), pid = replace_na(pid, 0),
         present = cov >= 0.5 & pid >= 50)

su <- full %>% group_by(pair_lab, cat) %>%
  summarise(n = n(), med_cov = median(cov), frac_present = mean(present), .groups = "drop")
w <- su %>% select(pair_lab, cat, frac_present) %>%
  pivot_wider(names_from = cat, values_from = frac_present)
cat("=== fraction of queries with an intact-looking match on the target ===\n")
print(as.data.frame(w), digits = 3)

cat("\n=== pooled by category ===\n")
print(as.data.frame(full %>% group_by(cat) %>%
  summarise(n = n(), med_cov = median(cov), frac_present = mean(present))), digits = 3)

if (all(c("lostX","lostY") %in% names(w))) {
  cat(sprintf("\nlostY > lostX in %d/8 | paired Wilcoxon p=%.3f\n",
              sum(w$lostY > w$lostX, na.rm = TRUE),
              suppressWarnings(wilcox.test(w$lostY, w$lostX, paired = TRUE)$p.value)))
  cat("lostY high = genes ARE on Y but unannotated => Y under-annotated, result in trouble\n")
  cat("lostY ~ lostX ~ negX = genes truly absent on both sides => real fractionation\n")
}
write_csv(full, "blast_presence_calls.csv"); write_csv(su, "blast_presence_summary.csv")
