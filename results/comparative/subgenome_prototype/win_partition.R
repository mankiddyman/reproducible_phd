#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Cluster each region's co-occurrence matrix into 2 groups = the subgenome split.
# Permutation null: shuffle tip->half labels within each window, preserving half sizes.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ape); library(ggplot2); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
set.seed(1); B <- 999

fs <- list.files("win/tre", pattern = "\\.treefile$", full.names = TRUE)
halves <- lapply(fs, function(f) {
  tr <- tryCatch(read.tree(f), error = function(e) NULL)
  if (is.null(tr) || !("NEP" %in% tr$tip.label) || Ntip(tr) < 5) return(NULL)
  tr <- root(tr, "NEP", resolve.root = TRUE); ing <- drop.tip(tr, "NEP")
  k <- ing$edge[ing$edge[, 1] == Ntip(ing) + 1, 2]
  if (length(k) != 2) return(NULL)
  h <- lapply(k, function(z) if (z <= Ntip(ing)) ing$tip.label[z] else extract.clade(ing, z)$tip.label)
  if (min(lengths(h)) < 1) return(NULL)
  list(region = sub("_w[0-9]+$", "", sub("\\.treefile$", "", basename(f))),
       win = sub("\\.treefile$", "", basename(f)), A = h[[1]], B = h[[2]])
})
halves <- Filter(Negate(is.null), halves)

pairmat <- function(hs, perm = FALSE) {
  tips <- sort(unique(unlist(lapply(hs, function(x) c(x$A, x$B)))))
  S <- N <- matrix(0, length(tips), length(tips), dimnames = list(tips, tips))
  for (x in hs) {
    all <- c(x$A, x$B)
    grp <- if (perm) { g <- sample(all); setNames(c(rep(1, length(x$A)), rep(2, length(x$B))), g) }
           else setNames(c(rep(1, length(x$A)), rep(2, length(x$B))), all)
    for (i in seq_along(all)) for (j in seq_along(all)) if (i < j) {
      a <- all[i]; b <- all[j]
      N[a, b] <- N[b, a] <- N[a, b] + 1
      if (grp[a] == grp[b]) S[a, b] <- S[b, a] <- S[a, b] + 1
    }
  }
  list(S = S, N = N, F = ifelse(N > 0, S / N, NA))
}

out <- list(); nullres <- list()
for (r in sort(unique(sapply(halves, `[[`, "region")))) {
  hs <- Filter(function(x) x$region == r, halves)
  if (length(hs) < 3) next
  m <- pairmat(hs)
  keep <- rownames(m$F)[rowSums(m$N >= 2, na.rm = TRUE) >= 2]
  if (length(keep) < 4) next
  F <- m$F[keep, keep]; F[is.na(F)] <- 0.5; diag(F) <- 1
  cl <- cutree(hclust(as.dist(1 - F), method = "average"), k = 2)
  # observed separation = mean within-group F minus mean between-group F
  sep <- function(F, cl) {
    w <- b <- c()
    for (i in seq_along(cl)) for (j in seq_along(cl)) if (i < j) {
      if (cl[i] == cl[j]) w <- c(w, F[i, j]) else b <- c(b, F[i, j]) }
    mean(w, na.rm = TRUE) - mean(b, na.rm = TRUE)
  }
  obs <- sep(F, cl)
  nul <- replicate(B, {
    mm <- pairmat(hs, perm = TRUE)
    FF <- mm$F[keep, keep]; FF[is.na(FF)] <- 0.5; diag(FF) <- 1
    sep(FF, cutree(hclust(as.dist(1 - FF), method = "average"), k = 2))
  })
  p <- (1 + sum(nul >= obs)) / (B + 1)
  g1 <- names(cl)[cl == 1]; g2 <- names(cl)[cl == 2]
  out[[r]] <- tibble(region = r, nwin = length(hs), ntax = length(keep),
                     separation = obs, null_mean = mean(nul), p = p,
                     group1 = paste(sort(g1), collapse = " "),
                     group2 = paste(sort(g2), collapse = " "))
}
res <- bind_rows(out) %>% mutate(p_adj = p.adjust(p, "BH"))
write_csv(res, "win_partitions.csv")

cat("=== inferred subgenome split per ancestral region ===\n")
for (i in seq_len(nrow(res))) {
  cat(sprintf("\n%-11s  windows=%d taxa=%d  separation=%.3f (null %.3f)  p=%.4f p_adj=%.4f\n",
              res$region[i], res$nwin[i], res$ntax[i], res$separation[i],
              res$null_mean[i], res$p[i], res$p_adj[i]))
  cat(sprintf("   A: %s\n   B: %s\n", res$group1[i], res$group2[i]))
}
cat(sprintf("\nregions with a significant split: %d/%d\n", sum(res$p_adj < 0.05), nrow(res)))

## does the split separate the two Dionaea homeologs, and which side is X?
ref <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE) %>%
  filter(anchor == "Nepenthes_gracilis")
xy <- bind_rows(transmute(ref, chr = winner, side = "X"),
                transmute(ref, chr = ifelse(winner == chrA, chrB, chrA), side = "Y")) %>%
  mutate(tag = paste0("Dio_", sub("^chr([0-9]+)_", "c\\1_", chr)))
chk <- res %>% rowwise() %>%
  mutate(dioA = paste(grep("^Dio", strsplit(group1, " ")[[1]], value = TRUE), collapse = ","),
         dioB = paste(grep("^Dio", strsplit(group2, " ")[[1]], value = TRUE), collapse = ","),
         sideA = paste(xy$side[match(strsplit(dioA, ",")[[1]], xy$tag)], collapse = ","),
         sideB = paste(xy$side[match(strsplit(dioB, ",")[[1]], xy$tag)], collapse = ",")) %>%
  ungroup() %>% select(region, dioA, sideA, dioB, sideB, p_adj)
cat("\n=== do the topological groups match the fractionation X/Y partition? ===\n")
print(as.data.frame(chk), digits = 3)
cat("one X on one side and one Y on the other in every region => the two signals agree\n")
