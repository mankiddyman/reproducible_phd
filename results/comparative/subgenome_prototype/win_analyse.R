#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# 60 window trees. Each window = one Nepenthes interval, all taxa share all genes.
# Q1: are the two Dionaea homeologs sisters, or interleaved with Drosera?
# Q2: does the partition REPEAT across consecutive windows of the same region?
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ape); library(ggplot2); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

fs <- list.files("win/tre", pattern = "\\.treefile$", full.names = TRUE)
cat(sprintf("window trees: %d\n", length(fs)))

info <- lapply(fs, function(f) {
  w <- sub("\\.treefile$", "", basename(f))
  tr <- tryCatch(read.tree(f), error = function(e) NULL)
  if (is.null(tr) || !("NEP" %in% tr$tip.label) || Ntip(tr) < 5) return(NULL)
  tr <- root(tr, "NEP", resolve.root = TRUE)
  ing <- drop.tip(tr, "NEP")
  kids <- ing$edge[ing$edge[, 1] == Ntip(ing) + 1, 2]
  if (length(kids) != 2) return(NULL)
  hv <- lapply(kids, function(k) if (k <= Ntip(ing)) ing$tip.label[k]
                                 else extract.clade(ing, k)$tip.label)
  dio <- grep("^Dio", tr$tip.label, value = TRUE)
  sup <- suppressWarnings(as.numeric(sub(".*/", "", tr$node.label)))
  # nearest heterospecific partner for each Dionaea copy
  nd <- node.depth(tr)
  part <- if (length(dio) >= 1) sapply(dio, function(d) {
    o <- setdiff(tr$tip.label, c(dio, "NEP"))
    if (!length(o)) return(NA_character_)
    o[which.min(sapply(o, function(x) nd[getMRCA(tr, c(d, x))]))]
  }) else character(0)
  list(win = w,
       region = sub("_w[0-9]+$", "", w),
       widx = as.integer(sub(".*_w", "", w)),
       ntip = Ntip(tr), ndio = length(dio),
       dio_sisters = if (length(dio) == 2) is.monophyletic(tr, dio) else NA,
       root_sup = if (all(is.na(sup))) NA_real_ else max(sup, na.rm = TRUE),
       med_sup = if (all(is.na(sup))) NA_real_ else median(sup, na.rm = TRUE),
       halves = hv, tips = tr$tip.label,
       partners = if (length(part)) paste(names(part), part, sep = "->", collapse = "; ") else NA)
})
info <- Filter(Negate(is.null), info)

sm <- bind_rows(lapply(info, function(x)
  tibble(win = x$win, region = x$region, widx = x$widx, ntip = x$ntip, ndio = x$ndio,
         dio_sisters = x$dio_sisters, root_sup = x$root_sup, med_sup = x$med_sup,
         halfA = paste(sort(x$halves[[1]]), collapse = ","),
         halfB = paste(sort(x$halves[[2]]), collapse = ","),
         partners = x$partners)))
write_csv(sm, "win_summary.csv")

cat("\n=== Q1: are the two Dionaea homeologs sisters? ===\n")
print(as.data.frame(sm %>% filter(ndio == 2) %>% group_by(region) %>%
  summarise(n = n(), sisters = sum(dio_sisters), frac = mean(dio_sisters),
            med_sup = round(median(med_sup, na.rm = TRUE)), .groups = "drop")))
q <- filter(sm, ndio == 2)
cat(sprintf("OVERALL: %d/%d windows have Dionaea copies as sisters (%.2f)\n",
            sum(q$dio_sisters), nrow(q), mean(q$dio_sisters)))
cat("low fraction => copies interleave with Drosera => shared WGD / phasable\n")

cat("\n=== deepest split per window (in genomic order) ===\n")
for (r in sort(unique(sm$region))) {
  cat(sprintf("\n-- %s --\n", r))
  d <- filter(sm, region == r) %>% arrange(widx)
  for (i in seq_len(nrow(d)))
    cat(sprintf("  w%02d sup=%-3s  [%s] || [%s]\n", d$widx[i],
                ifelse(is.na(d$root_sup[i]), "?", d$root_sup[i]), d$halfA[i], d$halfB[i]))
}

## ---- Q2: co-occurrence (orientation-free consistency) ----------------------
co <- bind_rows(lapply(info, function(x) {
  tp <- setdiff(x$tips, "NEP")
  if (length(tp) < 3) return(NULL)
  cb <- combn(sort(tp), 2)
  tibble(region = x$region, widx = x$widx,
         t1 = cb[1, ], t2 = cb[2, ],
         same = mapply(function(a, b)
           (a %in% x$halves[[1]] && b %in% x$halves[[1]]) ||
           (a %in% x$halves[[2]] && b %in% x$halves[[2]]), cb[1, ], cb[2, ]))
}))
cc <- co %>% group_by(region, t1, t2) %>%
  summarise(n = n(), same = mean(same), .groups = "drop") %>% filter(n >= 3)
write_csv(cc, "win_cooccurrence.csv")

cat("\n=== Q2: pairs seen in >=3 windows, fraction in the SAME half ===\n")
cat("near 1 or near 0 = reproducible partition | near 0.5 = no consistent structure\n")
for (r in sort(unique(cc$region))) {
  cat(sprintf("\n-- %s --\n", r))
  print(as.data.frame(filter(cc, region == r) %>% arrange(desc(same))), digits = 3)
}
cat(sprintf("\nOVERALL: %.0f%% of pairs are >=0.8 or <=0.2 (decisive), %.0f%% in 0.4-0.6 (noise)\n",
            100 * mean(cc$same >= 0.8 | cc$same <= 0.2), 100 * mean(cc$same > 0.4 & cc$same < 0.6)))

## ---- plot: consistency heatmap per region ---------------------------------
pl <- cc %>% mutate(pair = paste(t1, t2, sep = " | "))
p <- ggplot(pl, aes(same, reorder(pair, same), fill = same)) +
  geom_col() + geom_vline(xintercept = 0.5, linetype = 2) +
  facet_wrap(~ region, scales = "free_y", ncol = 2) +
  scale_fill_gradient2(low = "#C0392B", mid = "grey85", high = "#12795E", midpoint = 0.5) +
  labs(title = "Do two chromosomes land in the same half of the tree?",
       subtitle = "1 = always together, 0 = always apart, 0.5 = random. Across windows of one region.",
       x = "fraction of windows in the same half", y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "none")
ggsave("FIG9_window_consistency.png", p, width = 13, height = 14, dpi = 160, device = agg_png)
cat("\nWROTE: FIG9_window_consistency.png win_summary.csv win_cooccurrence.csv\n")
