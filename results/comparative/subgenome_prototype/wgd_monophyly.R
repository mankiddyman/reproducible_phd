#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Rooted on Nepenthes. For each species with >=2 copies: do its copies form a clade?
#   monophyletic => WGD postdates that species' divergence (independent)
#   interleaved  => WGD predates it (shared) -- and the pattern phases both species
# NULL: chance monophyly is NOT 0.5. Per-tree label permutation gives the real null.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(ape) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
B <- 200; set.seed(1)

meta <- read_tsv("wgd/tip_meta.tsv", show_col_types = FALSE)
fs <- list.files("wgd/tre", pattern = "\\.tre$", full.names = TRUE)
cat(sprintf("trees found: %d\n", length(fs)))

one <- function(f) {
  tr <- tryCatch(read.tree(f), error = function(e) NULL)
  if (is.null(tr) || Ntip(tr) < 4) return(NULL)
  m <- meta[match(tr$tip.label, meta$tip), ]
  if (any(is.na(m$genome))) return(NULL)
  nep <- tr$tip.label[m$genome == "Nepenthes_gracilis"]
  if (length(nep) != 1) return(NULL)
  tr <- tryCatch(root(tr, outgroup = nep, resolve.root = TRUE), error = function(e) NULL)
  if (is.null(tr)) return(NULL)
  m <- meta[match(tr$tip.label, meta$tip), ]
  a <- sub("\\.tre$", "", basename(f))
  n <- Ntip(tr)
  sup <- suppressWarnings(as.numeric(tr$node.label))
  msup <- if (all(is.na(sup))) NA_real_ else median(sup, na.rm = TRUE)
  g <- m$genome

  mono <- do.call(rbind, lapply(setdiff(unique(g), "Nepenthes_gracilis"), function(s) {
    idx <- which(g == s); k <- length(idx)
    if (k < 2) return(NULL)
    obs <- is.monophyletic(tr, tr$tip.label[idx])
    nullp <- mean(replicate(B, {
      p <- sample(setdiff(seq_len(n), which(g == "Nepenthes_gracilis")), k)
      is.monophyletic(tr, tr$tip.label[p])
    }))
    intr <- NA_character_
    if (!obs) {
      d <- extract.clade(tr, getMRCA(tr, tr$tip.label[idx]))
      o <- setdiff(unique(meta$genome[match(d$tip.label, meta$tip)]),
                   c(s, "Nepenthes_gracilis"))
      intr <- if (length(o)) paste(sort(o), collapse = "+") else "unresolved"
    }
    data.frame(anchor = a, sp = s, k = k, ntip = n, med_sup = msup,
               mono = obs, nullp = nullp, intruder = intr, stringsAsFactors = FALSE)
  }))

  nd <- node.depth(tr)
  ii <- which(g != "Nepenthes_gracilis")
  sis <- do.call(rbind, lapply(ii, function(i) {
    o <- ii[g[ii] != g[i]]
    if (!length(o)) return(NULL)
    dep <- vapply(o, function(j) nd[getMRCA(tr, tr$tip.label[c(i, j)])], 1)
    j <- o[which.min(dep)]
    data.frame(anchor = a, sp1 = g[i], chr1 = m$chr[i],
               sp2 = g[j], chr2 = m$chr[j], med_sup = msup, stringsAsFactors = FALSE)
  }))
  list(mono = mono, sis = sis)
}

r <- lapply(fs, one)
mono <- bind_rows(lapply(r, `[[`, "mono"))
sis  <- bind_rows(lapply(r, `[[`, "sis"))
write_csv(mono, "wgd/monophyly.csv"); write_csv(sis, "wgd/nearest_partner.csv")

report <- function(d, tag) {
  if (!nrow(d)) return(invisible())
  s <- d %>% group_by(sp) %>%
    summarise(n = n(), obs = sum(mono), exp = sum(nullp),
              frac_obs = mean(mono), frac_null = mean(nullp),
              z = (sum(mono) - sum(nullp)) / sqrt(sum(nullp * (1 - nullp))),
              .groups = "drop") %>%
    mutate(p = 2 * pnorm(-abs(z)))
  cat(sprintf("\n=== monophyly vs permutation null [%s] ===\n", tag))
  print(as.data.frame(s), digits = 3)
}
report(mono, "all trees")
report(filter(mono, !is.na(med_sup), med_sup >= 0.7), "median support >= 0.7")
cat("\nfrac_obs >> frac_null  => independent WGD\n")
cat("frac_obs ~ frac_null   => copies interleave => shared WGD\n")

cat("\n=== when not monophyletic, who nests inside? ===\n")
print(as.data.frame(mono %>% filter(!mono) %>% count(sp, intruder) %>%
  group_by(sp) %>% mutate(frac = round(n / sum(n), 3)) %>% ungroup()))

cat("\n=== nearest heterospecific partner ===\n")
print(as.data.frame(sis %>% count(sp1, sp2) %>% group_by(sp1) %>%
  mutate(frac = round(n / sum(n), 3)) %>% ungroup()))

dr <- sis %>% filter(sp1 == "Drosera_regia", sp2 == "Dionaea_muscipula") %>%
  count(chr1, chr2) %>% group_by(chr1) %>% mutate(tot = sum(n), frac = n / tot) %>%
  filter(tot >= 15) %>% slice_max(frac, n = 2) %>% ungroup() %>% arrange(chr1)
if (nrow(dr)) {
  cat("\n=== D. regia chromosomes: preferred Dionaea partner ===\n")
  print(as.data.frame(dr), digits = 3)
}
