#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Per ancestral region (Nepenthes chromosome): NJ tree of the chromosomes that
# descend from it. Distances are per-tree normalised so gene rate cancels.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ape); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
MINOBS <- 15

meta <- read_tsv("wgd/tip_meta.tsv", show_col_types = FALSE)
hits <- read_csv("synteny_ortho_hits.csv", show_col_types = FALSE)
anch <- hits %>% distinct(nep_gene, nep_chr)
ref  <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE) %>%
  filter(anchor == "Nepenthes_gracilis")
xy <- bind_rows(transmute(ref, chr = winner, side = "X"),
                transmute(ref, chr = ifelse(winner == chrA, chrB, chrA), side = "Y"))

lab <- meta %>% left_join(xy, by = "chr") %>%
  mutate(short = case_when(
    genome == "Nepenthes_gracilis" ~ "NEP",
    genome == "Dionaea_muscipula"  ~ paste0("Dio_", sub("^chr[0-9]+_", "", chr), "[", side, "]"),
    genome == "Drosera_regia"      ~ paste0("reg_", sub("_collapsed$", "", chr)),
    genome == "Drosera_binata"     ~ paste0("bin_", sub("_hap1$", "", chr)),
    TRUE ~ chr))

fs <- list.files("wgd/tre", pattern = "\\.tre$", full.names = TRUE)
acc <- new.env(hash = TRUE)

for (f in fs) {
  tr <- tryCatch(read.tree(f), error = function(e) NULL)
  if (is.null(tr) || Ntip(tr) < 4) next
  a <- sub("\\.tre$", "", basename(f))
  nc <- anch$nep_chr[match(a, anch$nep_gene)]
  if (is.na(nc)) next
  m <- lab[match(tr$tip.label, lab$tip), ]
  if (any(is.na(m$short))) next
  D <- cophenetic(tr)
  s <- mean(D[m$genome == "Nepenthes_gracilis", m$genome != "Nepenthes_gracilis"])
  if (!is.finite(s) || s <= 0) next
  D <- D / s
  for (i in seq_len(nrow(D) - 1)) for (j in (i + 1):nrow(D)) {
    ci <- m$short[i]; cj <- m$short[j]; if (ci == cj) next
    k <- paste(nc, min(ci, cj), max(ci, cj), sep = "|")
    v <- get0(k, envir = acc, ifnotfound = c(0, 0))
    assign(k, c(v[1] + D[i, j], v[2] + 1), envir = acc)
  }
}

df <- tibble(key = ls(acc)) %>%
  mutate(v = lapply(key, get, envir = acc),
         sum = vapply(v, `[`, 1, 1), n = vapply(v, `[`, 1, 2), d = sum / n) %>%
  separate(key, c("region", "c1", "c2"), sep = "\\|") %>% filter(n >= MINOBS)
write_csv(select(df, region, c1, c2, n, d), "region_chr_distances.csv")
cat(sprintf("chromosome pairs with >= %d observations: %d across %d regions\n",
            MINOBS, nrow(df), n_distinct(df$region)))

regs <- sort(unique(df$region))
agg_png("FIG3_region_trees.png", width = 1500, height = 850, res = 110)
par(mfrow = c(2, 4), mar = c(1, 1, 2.5, 1))
built <- 0
for (r in regs) {
  s <- filter(df, region == r)
  chrs <- sort(unique(c(s$c1, s$c2)))
  M <- matrix(NA_real_, length(chrs), length(chrs), dimnames = list(chrs, chrs))
  M[cbind(s$c1, s$c2)] <- s$d; M[cbind(s$c2, s$c1)] <- s$d; diag(M) <- 0
  if (anyNA(M)) {
    keep <- names(sort(colSums(is.na(M))))[1:max(3, sum(colSums(is.na(M)) == 0))]
    M <- M[keep, keep, drop = FALSE]
  }
  if (nrow(M) < 4 || anyNA(M)) { plot.new(); title(paste0(r, "\n(insufficient)")); next }
  tr <- tryCatch(nj(as.dist(M)), error = function(e) NULL)
  if (is.null(tr)) { plot.new(); title(r); next }
  if ("NEP" %in% tr$tip.label) tr <- root(tr, "NEP", resolve.root = TRUE)
  col <- ifelse(grepl("^Dio", tr$tip.label), "#d95f02",
         ifelse(grepl("^reg", tr$tip.label), "#7570b3",
         ifelse(grepl("^bin", tr$tip.label), "#1b9e77", "black")))
  plot(tr, tip.color = col, cex = 0.85, no.margin = FALSE, edge.width = 1.4)
  title(paste0(r, "  (n=", nrow(M), " chr)"), cex.main = 1)
  built <- built + 1
}
dev.off()
cat(sprintf("built %d region trees -> FIG3_region_trees.png\n", built))

## regia subgenome preference, restricted to informative trees
sis <- read_csv("wgd/nearest_partner.csv", show_col_types = FALSE)
ndio <- meta %>% filter(genome == "Dionaea_muscipula") %>% count(nep_gene, name = "kdio")
inf <- sis %>% filter(sp1 == "Drosera_regia", sp2 == "Dionaea_muscipula") %>%
  left_join(ndio, by = c("anchor" = "nep_gene")) %>% filter(kdio == 2) %>%
  left_join(xy, by = c("chr2" = "chr")) %>% filter(!is.na(side)) %>%
  count(chr1, side) %>% pivot_wider(names_from = side, values_from = n, values_fill = 0) %>%
  mutate(tot = X + Y, frac_X = X / tot) %>% filter(tot >= 20) %>%
  rowwise() %>% mutate(p = binom.test(X, tot, 0.5)$p.value) %>% ungroup() %>%
  mutate(p_adj = p.adjust(p, "BH")) %>% arrange(desc(frac_X))
cat("\n=== regia chromosomes: X vs Y preference (2-copy Dionaea trees only) ===\n")
print(as.data.frame(inf), digits = 3)
write_csv(inf, "regia_subgenome_preference.csv")
