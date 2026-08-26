#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Do complete-copy-number gene trees agree on a topology?
# Rooted on Nepenthes. Per tree: is each species' copy set monophyletic, and what
# is the CHROMOSOME composition of the two halves of the deepest ingroup split?
# Chromosome labels transfer between trees WITHIN an ancestral region, so a
# recurring chromosome partition is the subgenome signal.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ape); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

meta <- read_tsv("full14/tip_meta.tsv", show_col_types = FALSE)
if (!"nep_chr" %in% names(meta)) {
  cand <- grep("^nep_chr", names(meta), value = TRUE)
  if (length(cand)) meta$nep_chr <- meta[[cand[1]]]
}
if (!"nep_chr" %in% names(meta)) {
  h <- read_csv("synteny_ortho_hits.csv", show_col_types = FALSE) %>% distinct(nep_gene, nep_chr)
  meta <- left_join(meta, h, by = "nep_gene")
}
fullset <- readLines("full14/list_full.txt")
meta <- meta %>% mutate(
  set = ifelse(nep_gene %in% fullset, "full", "near"),
  lab = ifelse(genome == "Nepenthes_gracilis", "NEP",
        paste0(substr(sub("^Dro.era_|^Dionaea_", "", genome), 1, 3), "_",
               sub("_hap1$|_collapsed$", "", sub("^chr", "", chr)))))
cat(sprintf("meta: %d tips | full anchors %d | near anchors %d\n", nrow(meta),
            n_distinct(meta$nep_gene[meta$set == "full"]),
            n_distinct(meta$nep_gene[meta$set == "near"])))

fs <- list.files("full14/tre", pattern = "\\.treefile$", full.names = TRUE)

one <- function(f) {
  a <- sub("\\.treefile$", "", basename(f))
  tr <- tryCatch(read.tree(f), error = function(e) NULL)
  if (is.null(tr) || Ntip(tr) < 5) return(NULL)
  key <- function(x) gsub("@", "_", gsub("['\"]", "", x))
  tr$tip.label <- make.unique(tr$tip.label)
  m <- meta[match(key(tr$tip.label), key(meta$tip)), ]
  if (any(is.na(m$genome))) {
    message(sprintf("SKIP %s: %d/%d tips unmatched; e.g. %s", a,
                    sum(is.na(m$genome)), Ntip(tr),
                    paste(head(tr$tip.label[is.na(m$genome)], 2), collapse=" ")))
    return(NULL)
  }
  nep <- tr$tip.label[m$genome == "Nepenthes_gracilis"]
  if (length(nep) != 1) return(NULL)
  tr <- tryCatch(root(tr, nep, resolve.root = TRUE), error = function(e) NULL)
  if (is.null(tr)) return(NULL)
  m <- meta[match(key(tr$tip.label), key(meta$tip)), ]
  mono <- vapply(c("Dionaea_muscipula","Drosera_regia","Drosera_binata"), function(s) {
    tp <- tr$tip.label[m$genome == s]
    if (length(tp) < 2) NA else is.monophyletic(tr, tp) }, logical(1))
  ing <- drop.tip(tr, nep)
  if (Ntip(ing) < 4) return(NULL)
  kids <- ing$edge[ing$edge[, 1] == Ntip(ing) + 1, 2]
  if (length(kids) != 2) return(NULL)
  hv <- lapply(kids, function(k) if (k <= Ntip(ing)) ing$tip.label[k]
                                 else extract.clade(ing, k)$tip.label)
  lb <- function(x) paste(sort(meta$lab[match(key(x), key(meta$tip))]), collapse = ",")
  sh <- function(x) paste(sort(table(sub("_.*", "", meta$lab[match(key(x), key(meta$tip))]))), collapse = "/")
  o <- order(lengths(hv), vapply(hv, lb, ""))
  tibble(anchor = a, set = m$set[1], region = m$nep_chr[1], ntip = Ntip(tr),
         dio_mono = mono[1], reg_mono = mono[2], bin_mono = mono[3],
         halfA = lb(hv[[o[1]]]), halfB = lb(hv[[o[2]]]),
         shapeA = sh(hv[[o[1]]]), shapeB = sh(hv[[o[2]]]))
}

r <- bind_rows(lapply(fs, one))
write_csv(r, "full14_topologies.csv")
cat(sprintf("analysed %d trees (full %d, near %d)\n", nrow(r),
            sum(r$set == "full"), sum(r$set == "near")))

stopifnot(nrow(r) > 0)
fu <- filter(r, set == "full")
cat("\n================ FULL SET ================\n")
print(as.data.frame(select(fu, anchor, region, ntip, dio_mono, reg_mono, bin_mono)))
cat("\n--- monophyly rates ---\n")
print(as.data.frame(summarise(fu, n = n(),
  dio = mean(dio_mono, na.rm = TRUE), reg = mean(reg_mono, na.rm = TRUE),
  bin = mean(bin_mono, na.rm = TRUE))), digits = 3)
cat("\n--- deepest ingroup split ---\n")
for (i in seq_len(nrow(fu)))
  cat(sprintf("%-24s %-10s [%s] || [%s]\n", fu$anchor[i], fu$region[i], fu$halfA[i], fu$halfB[i]))

nr <- filter(r, set == "near")
cat("\n================ NEAR SET ================\n")
print(as.data.frame(summarise(nr, n = n(),
  dio = mean(dio_mono, na.rm = TRUE), reg = mean(reg_mono, na.rm = TRUE),
  bin = mean(bin_mono, na.rm = TRUE))), digits = 3)
cat("\n--- most frequent species compositions of the split ---\n")
print(as.data.frame(head(count(r, shapeA, shapeB, sort = TRUE), 10)))
cat("\n--- chromosome partitions recurring within a region (>=2 trees) ---\n")
print(as.data.frame(head(count(r, region, halfA, halfB, sort = TRUE) %>% filter(n >= 2), 30)))

## does a given chromosome pair co-occur in the same half more often than chance?
pairs <- bind_rows(lapply(seq_len(nrow(r)), function(i) {
  h1 <- strsplit(r$halfA[i], ",")[[1]]; h2 <- strsplit(r$halfB[i], ",")[[1]]
  g <- function(x, s) if (length(x) < 2) NULL else
    as_tibble(t(combn(sort(x), 2))) %>% setNames(c("c1","c2")) %>%
      mutate(region = r$region[i], same = TRUE)
  bind_rows(g(h1), g(h2),
    expand_grid(c1 = h1, c2 = h2) %>%
      mutate(a = pmin(c1, c2), b = pmax(c1, c2)) %>%
      transmute(c1 = a, c2 = b, region = r$region[i], same = FALSE))
})) %>% filter(c1 != c2, !grepl("^NEP", c1), !grepl("^NEP", c2))

co <- pairs %>% group_by(region, c1, c2) %>%
  summarise(n = n(), together = mean(same), .groups = "drop") %>%
  filter(n >= 4) %>% arrange(region, desc(together))
write_csv(co, "chromosome_cooccurrence.csv")
cat("\n--- chromosome pairs: how often in the SAME half (>=4 trees) ---\n")
print(as.data.frame(co), digits = 3)

agg_png("FIG7_full14_trees.png", width = 1700, height = 1050, res = 105)
par(mfrow = c(3, 5), mar = c(1, 1, 2.5, 1))
for (a in fu$anchor) {
  tr <- read.tree(file.path("full14/tre", paste0(a, ".treefile")))
  key <- function(x) gsub("@", "_", gsub("['\"]", "", x))
  m <- meta[match(key(tr$tip.label), key(meta$tip)), ]
  tr <- root(tr, tr$tip.label[m$genome == "Nepenthes_gracilis"], resolve.root = TRUE)
  m <- meta[match(key(tr$tip.label), key(meta$tip)), ]
  tr$tip.label <- m$lab
  col <- ifelse(grepl("^mus", tr$tip.label), "#d95f02",
         ifelse(grepl("^reg", tr$tip.label), "#7570b3",
         ifelse(grepl("^bin", tr$tip.label), "#1b9e77", "black")))
  plot(tr, tip.color = col, cex = 0.9, edge.width = 1.3)
  if (!is.null(tr$node.label))
    nodelabels(tr$node.label, frame = "none", cex = 0.6, adj = c(1.1, -0.4))
  title(sub("\\.t1$", "", a), cex.main = 0.9)
}
dev.off()
cat("\nWROTE: FIG7_full14_trees.png full14_topologies.csv chromosome_cooccurrence.csv\n")
