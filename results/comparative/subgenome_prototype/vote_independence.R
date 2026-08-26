#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# If the Drosera tips form a clade with both Dionaea copies outside it, every Drosera
# tip has the same MRCA with each Dionaea copy and they ALL vote identically by
# construction. Those trees carry one bit, not one bit per lineage.
# Measure how often that happens, then redo the key statistics on informative trees only.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(ape) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
set.seed(1); B <- 999

meta <- read_tsv("wgd7/tip_meta.tsv", show_col_types = FALSE)
key <- function(x) gsub("@", "_", gsub("['\"]", "", x))
v <- read_csv("tract_votes_blocks7.csv", show_col_types = FALSE) %>%
  mutate(species = sub("Drosera_","",sp))

## per tree: is the Drosera set monophyletic w.r.t. the two Dionaea tips?
res <- bind_rows(lapply(list.files("wgd7/tre", "\\.tre$", full.names = TRUE), function(f) {
  tr <- tryCatch(read.tree(f), error = function(e) NULL)
  if (is.null(tr) || Ntip(tr) < 4) return(NULL)
  m <- meta[match(key(tr$tip.label), key(meta$tip)), ]
  if (any(is.na(m$genome))) return(NULL)
  nep <- tr$tip.label[m$genome == "Nepenthes_gracilis"]
  if (length(nep) != 1) return(NULL)
  tr <- tryCatch(root(tr, nep, resolve.root = TRUE), error = function(e) NULL)
  if (is.null(tr)) return(NULL)
  m <- meta[match(key(tr$tip.label), key(meta$tip)), ]
  dio <- tr$tip.label[m$genome == "Dionaea_muscipula"]
  dro <- tr$tip.label[grepl("^Drosera_", m$genome)]
  if (length(dio) != 2 || length(dro) < 2) return(NULL)
  tibble(anchor = sub("\\.tre$", "", basename(f)),
         n_dro = length(dro), n_sp = n_distinct(m$genome[grepl("^Drosera_", m$genome)]),
         dro_clade = is.monophyletic(tr, dro))
}))
cat(sprintf("trees assessed: %d\n", nrow(res)))
cat(sprintf("Drosera tips form a CLADE (votes forced identical): %d (%.0f%%)\n",
            sum(res$dro_clade), 100*mean(res$dro_clade)))
cat("\nby number of Drosera species in the tree:\n")
print(as.data.frame(res %>% group_by(n_sp) %>%
  summarise(trees = n(), forced = sum(dro_clade), pct = round(100*mean(dro_clade)), .groups="drop")))

## do observed votes actually agree in the forced trees?
va <- v %>% inner_join(res, by = "anchor") %>%
  group_by(anchor, dro_clade) %>%
  summarise(n = n(), unanimous = n_distinct(vote) == 1, .groups="drop")
cat("\n=== unanimity of votes, by whether Drosera is a clade ===\n")
print(as.data.frame(va %>% group_by(dro_clade) %>%
  summarise(anchors = n(), unanimous = sum(unanimous),
            pct = round(100*mean(unanimous)), .groups="drop")))

## ---- redo pairwise agreement on INFORMATIVE trees only --------------------
SPO <- c("regia","binata","paradoxa","scorpioides","capensis")
agree_tab <- function(d, label) {
  call <- d %>% group_by(pair, anchor, species) %>%
    summarise(k = n(), a = sum(vote == "A"), .groups="drop") %>%
    filter(a != k - a) %>% mutate(call = ifelse(a > k - a, 1L, 0L))
  W <- call %>% select(pair, anchor, species, call) %>%
    pivot_wider(names_from = species, values_from = call)
  sp <- intersect(SPO, names(W))
  m <- e <- matrix(NA_real_, length(sp), length(sp), dimnames = list(sp, sp))
  for (i in seq_along(sp)) for (j in seq_along(sp)) {
    x <- W[[sp[i]]]; y <- W[[sp[j]]]; ok <- !is.na(x) & !is.na(y)
    if (sum(ok) >= 20) { a <- mean(x[ok]==y[ok]); px <- mean(x[ok]); py <- mean(y[ok])
      m[i,j] <- a; e[i,j] <- a - (px*py + (1-px)*(1-py)) } }
  cat(sprintf("\n=== %s (anchors: %d) ===\nagreement:\n", label, n_distinct(call$anchor)))
  print(round(m,3)); cat("excess:\n"); print(round(e,3))
  invisible(list(agr=m, exc=e))
}
vi <- v %>% inner_join(filter(res, !dro_clade) %>% select(anchor), by = "anchor")
agree_tab(v,  "ALL trees")
agree_tab(vi, "INFORMATIVE trees only (Drosera not a clade)")

## ---- redo the concordance statistic on informative trees ------------------
conc <- function(d, label) {
  lin <- d %>% group_by(pair, anchor, sp, lineage) %>%
    summarise(fa = mean(vote=="A"), .groups="drop") %>%
    filter(fa != .5) %>% mutate(c = ifelse(fa > .5, 1L, 0L))
  g <- lin %>% group_by(pair, anchor) %>%
    summarise(k = n(), x = sum(c), .groups="drop") %>% filter(k >= 2) %>%
    mutate(agree = pmax(x, k-x)/k)
  out <- do.call(rbind, lapply(split(g, g$pair), function(dd) {
    pA <- sum(dd$x)/sum(dd$k)
    nul <- vapply(seq_len(B), function(i) { xx <- rbinom(nrow(dd), dd$k, pA)
      mean(pmax(xx, dd$k-xx)/dd$k) }, 0)
    data.frame(pair=dd$pair[1], n=nrow(dd), obs=mean(dd$agree), null=mean(nul),
               excess=mean(dd$agree)-mean(nul),
               z=(mean(dd$agree)-mean(nul))/sd(nul)) }))
  cat(sprintf("\n=== concordance, %s ===\n", label)); print(out, digits=3, row.names=FALSE)
  cat(sprintf("mean excess %.3f | mean z %.2f\n", mean(out$excess), mean(out$z)))
}
conc(v,  "ALL trees")
conc(vi, "INFORMATIVE trees only")
