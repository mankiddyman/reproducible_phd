#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Orient every window by the Dionaea X/Y pair (X/Y from fractionation = independent).
# Then every Drosera block votes X-side or Y-side. Pooling across all windows gives
# per-block power that per-region clustering cannot reach.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ape); library(ggplot2); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
set.seed(1); B <- 999

ref <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE) %>%
  filter(anchor == "Nepenthes_gracilis")
xy <- bind_rows(transmute(ref, chr = winner, side = "X"),
                transmute(ref, chr = ifelse(winner == chrA, chrB, chrA), side = "Y")) %>%
  mutate(tag = paste0("Dio_", sub("^chr([0-9]+)_", "c\\1_", chr)))

fs <- list.files("win25/tre", pattern = "\\.treefile$", full.names = TRUE)

W <- lapply(fs, function(f) {
  tr <- tryCatch(read.tree(f), error = function(e) NULL)
  if (is.null(tr) || !("NEP" %in% tr$tip.label) || Ntip(tr) < 5) return(NULL)
  tr <- root(tr, "NEP", resolve.root = TRUE); ing <- drop.tip(tr, "NEP")
  k <- ing$edge[ing$edge[, 1] == Ntip(ing) + 1, 2]
  if (length(k) != 2) return(NULL)
  h <- lapply(k, function(z) if (z <= Ntip(ing)) ing$tip.label[z] else extract.clade(ing, z)$tip.label)
  w <- sub("\\.treefile$", "", basename(f))
  dio <- grep("^Dio", ing$tip.label, value = TRUE)
  sx <- xy$tag[xy$side == "X"]; sy <- xy$tag[xy$side == "Y"]
  list(win = w, region = sub("_w[0-9]+$", "", w), tips = ing$tip.label, halves = h,
       dioX = intersect(dio, sx), dioY = intersect(dio, sy))
})
W <- Filter(Negate(is.null), W)
cat(sprintf("windows: %d | with both Dio copies: %d\n", length(W),
            sum(sapply(W, function(x) length(x$dioX) == 1 && length(x$dioY) == 1))))

## ---- TEST 1: are the Dionaea homeologs separated by the deepest split? -----
sep <- bind_rows(lapply(W, function(x) {
  if (length(x$dioX) != 1 || length(x$dioY) != 1) return(NULL)
  a <- length(x$halves[[1]]); b <- length(x$halves[[2]]); n <- a + b
  inA <- function(t) t %in% x$halves[[1]]
  tibble(win = x$win, region = x$region, ntip = n,
         separated = inA(x$dioX) != inA(x$dioY),
         exp = 2 * a * b / (n * (n - 1)))
}))
o <- sum(sep$separated); e <- sum(sep$exp)
z <- (o - e) / sqrt(sum(sep$exp * (1 - sep$exp)))
cat(sprintf("\n=== TEST 1: Dionaea homeologs on opposite sides ===\nobserved %d/%d, expected %.1f, z=%.2f, p=%.4g\n",
            o, nrow(sep), e, z, 2 * pnorm(-abs(z))))
print(as.data.frame(sep %>% group_by(region) %>%
  summarise(n = n(), sep = sum(separated), exp = round(sum(exp), 1), .groups = "drop")))

## ---- TEST 2: which side does each Drosera block sit on? -------------------
vot <- bind_rows(lapply(W, function(x) {
  if (length(x$dioX) != 1 || length(x$dioY) != 1) return(NULL)
  inA <- function(t) t %in% x$halves[[1]]
  if (inA(x$dioX) == inA(x$dioY)) return(NULL)     # unoriented window
  xa <- inA(x$dioX)
  t <- setdiff(x$tips, c(x$dioX, x$dioY))
  if (!length(t)) return(NULL)
  tibble(win = x$win, region = x$region, taxon = t,
         side = ifelse(inA(t) == xa, "X", "Y"))
}))
cat(sprintf("\noriented windows: %d | block votes: %d\n", n_distinct(vot$win), nrow(vot)))

asg <- vot %>% group_by(region, taxon) %>%
  summarise(n = n(), X = sum(side == "X"), fracX = mean(side == "X"), .groups = "drop") %>%
  filter(n >= 3) %>% rowwise() %>%
  mutate(p = binom.test(X, n, 0.5)$p.value) %>% ungroup() %>%
  mutate(p_adj = p.adjust(p, "BH"),
         call = case_when(p_adj < 0.1 & fracX > .5 ~ "X",
                          p_adj < 0.1 & fracX < .5 ~ "Y", TRUE ~ "?")) %>%
  arrange(region, desc(fracX))
write_csv(asg, "block_subgenome_calls_w25.csv")
cat("\n=== TEST 2: Drosera block -> subgenome ===\n")
print(as.data.frame(asg), digits = 3)
cat(sprintf("\ncalled: %d X, %d Y, %d ambiguous (of %d blocks)\n",
            sum(asg$call == "X"), sum(asg$call == "Y"), sum(asg$call == "?"), nrow(asg)))

## ---- TEST 3: is the whole thing better than chance? -----------------------
obs <- mean(abs(asg$fracX - 0.5))
nul <- replicate(B, {
  v <- vot %>% group_by(win) %>% mutate(side = sample(side)) %>% ungroup() %>%
    group_by(region, taxon) %>% summarise(n = n(), f = mean(side == "X"), .groups = "drop") %>%
    filter(n >= 3)
  mean(abs(v$f - 0.5))
})
cat(sprintf("\n=== TEST 3: global consistency ===\nmean |fracX-0.5| = %.3f | null %.3f (sd %.3f) | p = %.4f\n",
            obs, mean(nul), sd(nul), (1 + sum(nul >= obs)) / (B + 1)))

p <- ggplot(asg, aes(fracX, reorder(paste(region, taxon), fracX), fill = call)) +
  geom_col() + geom_vline(xintercept = .5, linetype = 2) +
  scale_fill_manual(values = c(X = "#12795E", Y = "#C0392B", `?` = "grey75")) +
  labs(title = "Drosera block assignment to the Dionaea subgenomes",
       subtitle = "windows oriented by the Dionaea X/Y pair; 1 = always with X, 0 = always with Y",
       x = "fraction of windows on the X side", y = NULL, fill = "call") +
  theme_bw(base_size = 9)
ggsave("FIG10_block_calls_w25.png", p, width = 11, height = 10, dpi = 170, device = agg_png)
cat("\nWROTE: FIG10_block_calls_w25.png block_subgenome_calls_w25.csv\n")
