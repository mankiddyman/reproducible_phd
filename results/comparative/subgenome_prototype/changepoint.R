#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Changepoint scan with permutation null.
# For each (Dionaea chromosome pair, Drosera lineage-block): order the votes along the
# Dionaea coordinate, scan every cut, score by binomial log-likelihood ratio, keep the
# best. Then shuffle the vote ORDER (same votes, no positional structure) and rescan,
# B times. p = fraction of shuffles whose best cut scores at least as high.
# The LLR is used rather than |pL - pR| because it penalises tiny segments.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ragg); library(patchwork)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
MINN <- 15; MINSEG <- 5; B <- 999
set.seed(1)

v <- read_csv("tract_votes_blocks.csv", show_col_types = FALSE) %>%
  mutate(unit = sub("Drosera_[a-z]+_vs_Nepenthes_gracilis: ", "", unit))

ll <- function(x, n) {
  if (n <= 0) return(0)
  r <- 0
  if (x > 0)  r <- r + x * log(x / n)
  if (x < n)  r <- r + (n - x) * log((n - x) / n)
  r
}
scan_llr <- function(b) {                      # b = 0/1 vector in genomic order
  n <- length(b); cs <- cumsum(b); tot <- cs[n]
  base <- ll(tot, n); ks <- MINSEG:(n - MINSEG)
  if (!length(ks)) return(c(llr = 0, k = NA_real_))
  val <- vapply(ks, function(k) ll(cs[k], k) + ll(tot - cs[k], n - k) - base, 0)
  i <- which.max(val); c(llr = val[i], k = ks[i])
}

S <- v %>% arrange(pair, unit, posA) %>%
  group_by(pair, chrA, chrB, unit) %>% filter(n() >= MINN) %>%
  summarise(n = n(), fracA = mean(vote == "A"),
            pos = list(posA), b = list(as.integer(vote == "A")), .groups = "drop")
cat(sprintf("series with >=%d votes: %d (across %d Dionaea pairs)\n",
            MINN, nrow(S), n_distinct(S$pair)))
if (!nrow(S)) { cat("nothing to scan\n"); quit(save = "no") }

res <- bind_rows(lapply(seq_len(nrow(S)), function(i) {
  b <- S$b[[i]]; pos <- S$pos[[i]]; n <- length(b)
  o <- scan_llr(b)
  if (is.na(o["k"])) return(NULL)
  nul <- replicate(B, scan_llr(sample(b))["llr"])
  k <- as.integer(o["k"])
  tibble(pair = S$pair[i], chrA = S$chrA[i], chrB = S$chrB[i], unit = S$unit[i],
         n = n, fracA = S$fracA[i], k = k,
         cut_Mb = (pos[k] + pos[k+1]) / 2 / 1e6,
         left_fracA  = mean(b[1:k]),
         right_fracA = mean(b[(k+1):n]),
         llr = unname(o["llr"]), null_q95 = unname(quantile(nul, .95)),
         p = (1 + sum(nul >= o["llr"])) / (B + 1))
}))
res <- res %>% mutate(p_adj = p.adjust(p, "BH"),
                      verdict = ifelse(p_adj < 0.05, "BREAK", "no break")) %>%
  arrange(p)
write_csv(res, "changepoint_results.csv")

cat("\n=== top series by evidence for a breakpoint ===\n")
print(as.data.frame(res %>% transmute(pair, unit, n, fracA = round(fracA,2),
  cut_Mb = round(cut_Mb,1), left = round(left_fracA,2), right = round(right_fracA,2),
  llr = round(llr,2), null95 = round(null_q95,2), p = round(p,4),
  p_adj = round(p_adj,3), verdict) %>% head(20)), right = FALSE)

cat(sprintf("\nsignificant breakpoints (BH<0.05): %d of %d series\n",
            sum(res$verdict == "BREAK"), nrow(res)))
cat(sprintf("p<0.05 uncorrected: %d observed vs %.1f expected by chance\n",
            sum(res$p < 0.05), 0.05 * nrow(res)))
cat(sprintf("median observed LLR %.2f vs median null 95th pct %.2f\n",
            median(res$llr), median(res$null_q95)))

## ---- plots ----------------------------------------------------------------
tracks <- bind_rows(lapply(seq_len(nrow(S)), function(i)
  tibble(pair = S$pair[i], unit = S$unit[i],
         pos_Mb = S$pos[[i]] / 1e6, vote = ifelse(S$b[[i]] == 1, "A", "B")))) %>%
  left_join(select(res, pair, unit, cut_Mb, verdict), by = c("pair", "unit")) %>%
  mutate(lab = paste(pair, unit))

p1 <- ggplot(tracks, aes(pos_Mb, lab, colour = vote)) +
  geom_point(shape = 124, size = 3, alpha = .85) +
  geom_point(data = distinct(tracks, pair, lab, cut_Mb, verdict) %>%
               filter(verdict == "BREAK"),
             aes(cut_Mb, lab), inherit.aes = FALSE, shape = 4, size = 3.5,
             stroke = 1.1, colour = "black") +
  facet_wrap(~ pair, scales = "free", ncol = 2) +
  scale_colour_manual(values = c(A = "#1b9e77", B = "#d95f02")) +
  labs(title = "Vote tracks along the Dionaea chromosome, one row per Drosera block",
       subtitle = "green = allies sgA, orange = allies sgB; black cross = significant changepoint (BH<0.05)",
       x = "position on sgA (Mb)", y = NULL, colour = "vote") +
  theme_bw(base_size = 8)
ggsave("FIG19_vote_tracks.png", p1, width = 13, height = 11, dpi = 170, device = agg_png)

p2 <- ggplot(res, aes(null_q95, llr, colour = verdict)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  geom_point(aes(size = n), alpha = .8) +
  scale_colour_manual(values = c(BREAK = "#d95f02", `no break` = "grey60")) +
  scale_size_continuous(range = c(1.5, 4)) +
  labs(title = "Observed best-cut score vs its own permutation null",
       subtitle = "dashed line = the 95th percentile of shuffles; points above it beat their null",
       x = "null 95th percentile LLR", y = "observed best-cut LLR",
       colour = NULL, size = "votes") +
  theme_bw(base_size = 10)

p3 <- ggplot(res, aes(p)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "#7570b3", colour = "white") +
  geom_hline(yintercept = nrow(res) * 0.05, linetype = 2, colour = "grey40") +
  labs(title = "Distribution of permutation p-values",
       subtitle = "dashed line = what a uniform (no-signal) distribution would give per bin",
       x = "permutation p", y = "series") +
  theme_bw(base_size = 10)

ggsave("FIG20_changepoint_null.png", p2 / p3, width = 8, height = 9, dpi = 180, device = agg_png)
cat("\nWROTE: FIG19_vote_tracks.png FIG20_changepoint_null.png changepoint_results.csv\n")
