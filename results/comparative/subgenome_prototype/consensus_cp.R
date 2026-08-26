#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# CONSENSUS across Drosera lineages, one track per Dionaea chromosome pair.
# Rationale: lineage-specific gene conversion / partial HE hits different loci in
# different Drosera lineages, so it averages out. A Dionaea HE flips EVERY lineage
# at the same coordinate, so it survives averaging.
# At each Dionaea anchor position: each Drosera lineage-block casts ONE vote
# (majority of its genes there; ties dropped). x = lineages voting sgA, k = lineages voting.
# Changepoint = binomial LLR on the (x, k) series, permutation null on anchor order.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ragg); library(patchwork)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
MINSEG <- 5; MINANCH <- 20; B <- 999
set.seed(1)

v <- read_csv("tract_votes_blocks6.csv", show_col_types = FALSE) %>%
  mutate(unit = sub("Drosera_[a-z]+_vs_Nepenthes_gracilis: ", "", unit))
cat(sprintf("raw votes: %d | anchors: %d | lineage-blocks: %d | Dionaea pairs: %d\n",
            nrow(v), n_distinct(v$anchor), n_distinct(v$unit), n_distinct(v$pair)))

## one vote per lineage-block per anchor
LV <- v %>% group_by(pair, chrA, chrB, anchor, posA, sp, unit) %>%
  summarise(fa = mean(vote == "A"), .groups = "drop") %>%
  filter(fa != 0.5) %>% mutate(call = ifelse(fa > 0.5, 1L, 0L))

## consensus per anchor
CS <- LV %>% group_by(pair, chrA, chrB, anchor, posA) %>%
  summarise(k = n(), x = sum(call), nsp = n_distinct(sp), .groups = "drop") %>%
  arrange(pair, posA)
cat(sprintf("\nconsensus positions: %d (vs %d series in the per-lineage version)\n", nrow(CS), 64))
cat("\nlineage-blocks voting per position:\n"); print(as.data.frame(count(CS, k)))
cat("\nspecies voting per position:\n"); print(as.data.frame(count(CS, nsp)))
cat("\npositions per Dionaea pair:\n")
print(as.data.frame(CS %>% group_by(pair) %>%
  summarise(n_pos = n(), median_k = median(k), total_votes = sum(k),
            span_Mb = round((max(posA)-min(posA))/1e6, 1), .groups = "drop")))

ll <- function(x, n) { if (n <= 0) return(0); r <- 0
  if (x > 0) r <- r + x*log(x/n); if (x < n) r <- r + (n-x)*log((n-x)/n); r }
scan_cnt <- function(x, k) {
  n <- length(x); cx <- cumsum(x); ck <- cumsum(k)
  base <- ll(cx[n], ck[n]); ks <- MINSEG:(n - MINSEG)
  if (!length(ks)) return(c(llr = 0, i = NA_real_))
  val <- vapply(ks, function(j)
    ll(cx[j], ck[j]) + ll(cx[n]-cx[j], ck[n]-ck[j]) - base, 0)
  m <- which.max(val); c(llr = val[m], i = ks[m])
}

S <- CS %>% group_by(pair, chrA, chrB) %>% filter(n() >= MINANCH) %>%
  summarise(n = n(), pos = list(posA), x = list(x), k = list(k), .groups = "drop")
cat(sprintf("\npairs with >=%d consensus positions: %d\n", MINANCH, nrow(S)))

res <- bind_rows(lapply(seq_len(nrow(S)), function(i) {
  x <- S$x[[i]]; k <- S$k[[i]]; pos <- S$pos[[i]]; n <- length(x)
  o <- scan_cnt(x, k); if (is.na(o["i"])) return(NULL)
  nul <- replicate(B, { p <- sample(n); scan_cnt(x[p], k[p])["llr"] })
  j <- as.integer(o["i"])
  tibble(pair = S$pair[i], chrA = S$chrA[i], chrB = S$chrB[i], n_pos = n,
         total_votes = sum(k), fracA = sum(x)/sum(k),
         cut_Mb = (pos[j] + pos[j+1])/2/1e6,
         left_fracA  = sum(x[1:j])/sum(k[1:j]),
         right_fracA = sum(x[(j+1):n])/sum(k[(j+1):n]),
         llr = unname(o["llr"]), null_q95 = unname(quantile(nul, .95)),
         p = (1 + sum(nul >= o["llr"]))/(B + 1))
})) %>% mutate(p_adj = p.adjust(p, "BH"),
               verdict = ifelse(p_adj < 0.05, "BREAK", "no break")) %>% arrange(p)
write_csv(res, "consensus_changepoints.csv")

cat("\n=== consensus changepoint per Dionaea chromosome pair ===\n")
print(as.data.frame(res %>% transmute(pair, n_pos, total_votes,
  fracA = round(fracA,3), cut_Mb = round(cut_Mb,1),
  left = round(left_fracA,2), right = round(right_fracA,2),
  llr = round(llr,2), null95 = round(null_q95,2),
  p = round(p,4), p_adj = round(p_adj,3), verdict)), right = FALSE)
cat(sprintf("\nsignificant: %d of %d pairs\n", sum(res$verdict == "BREAK"), nrow(res)))

## ---- plots ----------------------------------------------------------------
pl <- CS %>% inner_join(select(res, pair, cut_Mb, verdict), by = "pair") %>%
  mutate(fracA = x/k, pos_Mb = posA/1e6)
sm <- pl %>% group_by(pair) %>% arrange(pos_Mb) %>%
  mutate(roll = zoo_roll <- {
    w <- 9; n <- n()
    if (n < w) rep(NA_real_, n) else
      sapply(seq_len(n), function(i) { lo <- max(1, i-(w-1)%/%2); hi <- min(n, i+(w-1)%/%2)
        sum(x[lo:hi])/sum(k[lo:hi]) }) }) %>% ungroup()

p1 <- ggplot(sm, aes(pos_Mb, fracA)) +
  geom_hline(yintercept = .5, linetype = 2, colour = "grey40") +
  geom_point(aes(size = k), alpha = .45, colour = "#7570b3") +
  geom_line(aes(y = roll), colour = "#d95f02", linewidth = .8, na.rm = TRUE) +
  geom_vline(data = filter(res, verdict == "BREAK"), aes(xintercept = cut_Mb),
             colour = "black", linetype = 3) +
  facet_wrap(~ pair, ncol = 2, scales = "free_x") +
  scale_size_continuous(range = c(.8, 3.5)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Consensus across ALL Drosera lineages, per Dionaea chromosome pair",
       subtitle = "point = one Dionaea gene, size = lineages voting; orange = 9-position rolling consensus",
       x = "position on sgA (Mb)", y = "fraction of lineages voting sgA", size = "lineages") +
  theme_bw(base_size = 9)
ggsave("FIG21_consensus_tracks.png", p1, width = 12, height = 10, dpi = 175, device = agg_png)

p2 <- ggplot(res, aes(null_q95, llr, colour = verdict)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  geom_point(aes(size = total_votes)) +
  ggrepel::geom_text_repel(aes(label = pair), size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c(BREAK = "#d95f02", `no break` = "grey55")) +
  labs(title = "Consensus best-cut score vs its permutation null",
       x = "null 95th percentile LLR", y = "observed LLR", colour = NULL, size = "votes") +
  theme_bw(base_size = 10)
ggsave("FIG22_consensus_null.png", p2, width = 8, height = 5.5, dpi = 180, device = agg_png)
cat("\nWROTE: FIG21_consensus_tracks.png FIG22_consensus_null.png consensus_changepoints.csv\n")
