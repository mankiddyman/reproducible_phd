#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Per-lineage skew is the unit. Each Drosera chromosome is scored independently
# against the two Dionaea homeologs, so unequal retention/conversion between
# lineages is irrelevant. Then propagate labels across Dionaea pairs via shared
# Drosera chromosomes and test whether the constraint graph two-colours.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
MINN <- 15; ALPHA <- 0.05

v <- read_csv("tract_votes.csv", show_col_types = FALSE)

## ---- 1. per-lineage skew ----------------------------------------------------
L <- v %>% group_by(pair, chrA, chrB, lineage) %>%
  summarise(n = n(), A = sum(vote == "A"), frac_A = A/n,
            span_Mb = (max(posA) - min(posA))/1e6, .groups = "drop") %>%
  filter(n >= MINN) %>% rowwise() %>%
  mutate(p = binom.test(A, n, 0.5)$p.value,
         lo = binom.test(A, n, 0.5)$conf.int[1],
         hi = binom.test(A, n, 0.5)$conf.int[2]) %>% ungroup() %>%
  group_by(pair) %>% mutate(p_adj = p.adjust(p, "BH")) %>% ungroup() %>%
  mutate(call = case_when(p_adj < ALPHA & frac_A > .5 ~ "A",
                          p_adj < ALPHA & frac_A < .5 ~ "B", TRUE ~ "?")) %>%
  arrange(pair, desc(frac_A))
write_csv(L, "lineage_skew.csv")

cat("=== per-lineage skew, by Dionaea pair ===\n")
for (p in unique(L$pair)) {
  d <- filter(L, pair == p)
  cat(sprintf("\n-- %s   (%s = A / %s = B) --\n", p, d$chrA[1], d$chrB[1]))
  print(as.data.frame(select(d, lineage, n, frac_A, lo, hi, p_adj, call)), digits = 3)
}
cat(sprintf("\ncalled: %d A, %d B, %d ambiguous (of %d lineages)\n",
            sum(L$call=="A"), sum(L$call=="B"), sum(L$call=="?"), nrow(L)))

## is the split within a pair genuinely bimodal?
cat("\n=== bimodality within each pair ===\n")
print(as.data.frame(L %>% group_by(pair) %>%
  summarise(n_lin = n(), nA = sum(call=="A"), nB = sum(call=="B"),
            spread = round(max(frac_A) - min(frac_A), 3),
            both_sides = nA > 0 & nB > 0, .groups = "drop")))

## ---- 2. cross-pair constraint graph ----------------------------------------
called <- filter(L, call != "?") %>%
  mutate(sp_chr = lineage)   # a Drosera chromosome, shared across Dionaea pairs

links <- called %>% group_by(sp_chr) %>% filter(n_distinct(pair) >= 2) %>%
  arrange(pair) %>%
  summarise(pairs = list(pair), calls = list(call), n = n(), .groups = "drop")
cat(sprintf("\n=== Drosera chromosomes linking >=2 Dionaea pairs: %d ===\n", nrow(links)))
if (nrow(links)) {
  ed <- bind_rows(lapply(seq_len(nrow(links)), function(i) {
    ps <- links$pairs[[i]]; cs <- links$calls[[i]]
    cb <- combn(seq_along(ps), 2)
    tibble(sp_chr = links$sp_chr[i], p1 = ps[cb[1,]], p2 = ps[cb[2,]],
           same = cs[cb[1,]] == cs[cb[2,]])
  }))
  print(as.data.frame(ed))

  ## greedy two-colouring, count violated constraints
  nodes <- sort(unique(c(ed$p1, ed$p2)))
  col <- setNames(rep(NA_integer_, length(nodes)), nodes)
  col[nodes[1]] <- 1L
  for (it in 1:20) for (i in seq_len(nrow(ed))) {
    a <- ed$p1[i]; b <- ed$p2[i]; s <- ed$same[i]
    if (!is.na(col[a]) && is.na(col[b])) col[b] <- if (s) col[a] else 3L - col[a]
    if (!is.na(col[b]) && is.na(col[a])) col[a] <- if (s) col[b] else 3L - col[b]
  }
  viol <- sum(mapply(function(a,b,s) {
    if (is.na(col[a]) || is.na(col[b])) return(FALSE)
    (col[a] == col[b]) != s }, ed$p1, ed$p2, ed$same))
  cat(sprintf("\n=== two-colouring: %d/%d constraints violated ===\n", viol, nrow(ed)))
  print(data.frame(pair = names(col), group = ifelse(col == 1, "A", "B")))
  cat("0 violations => a single global partition is consistent with the links\n")
  cat("many violations => links unreliable (mosaic Drosera chromosomes) or no global partition\n")
} else cat("no shared lineages -> cannot link pairs; global orientation not derivable this way\n")

## ---- 3. plot ---------------------------------------------------------------
p1 <- ggplot(L, aes(frac_A, reorder(paste(pair, lineage), frac_A), colour = call)) +
  geom_vline(xintercept = .5, linetype = 2, colour = "grey40") +
  geom_linerange(aes(xmin = lo, xmax = hi)) +
  geom_point(aes(size = n)) +
  facet_wrap(~ pair, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c(A = "#1b9e77", B = "#d95f02", `?` = "grey65")) +
  scale_size_continuous(range = c(1.2, 3.2)) +
  labs(title = "Which Dionaea homeolog does each Drosera chromosome ally with?",
       subtitle = "one row per Drosera chromosome; 95% CI; each lineage scored independently",
       x = "fraction of genes voting sgA", y = NULL, colour = "call", size = "genes") +
  theme_bw(base_size = 9)
ggsave("FIG14_lineage_skew.png", p1, width = 11, height = 10, dpi = 170, device = agg_png)
cat("\nWROTE: FIG14_lineage_skew.png lineage_skew.csv\n")
