#!/usr/bin/env Rscript
# =============================================================================
# DR21c - ARE THE COPIES AT ONE LOCUS LABELLED TOGETHER OR INDEPENDENTLY?
# =============================================================================
# A hexaploid locus that kept all three copies should read exactly A, A, B under
# AAB. If labels were independent draws at the species' A-rate pA, the three
# would read AAA / AAB / ABB / BBB binomially. AAB far above chance and the rest
# far below means the copies at a locus share a fixed composition: the
# constitution.
#
# What this adds to the original pooled figure
#   PER SPECIES   the pooled figure mixes five species, capensis included, whose
#                 3-copy loci are not full retention (3 of its 6 copies).
#   HONEST NULL   labels are set per SEGMENT (DR02: sign of the segment's median
#                 delta), so every 3-copy locus in a region reuses the same few
#                 segments and the loci are not independent replicates. Null:
#                 shuffle segment labels within a species (keeping how many are
#                 A and B), recompute every locus, 2000 times.
#   2-COPY LOCI   under AAB a 2-copy locus is AA or AB, never BB (script 29's
#                 "A = 0" test), on several times more loci.
#   BY REGION     a Dionaea pair with the wrong A/B orientation would make that
#                 region read ABB in every species.
#   BUBBLE LINK   the 3-copy figure is the A + B = 3 diagonal of DR21's bubbles.
#
# Inputs  DR/out/DR21_copynumber.csv, DR02_propagated_blocks.csv, pairwise_ks.csv
# Outputs DR/out/DR21c_*.csv, DR/fig/DR21c_*.png, DR/fig/DR21c_all.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})
options(dplyr.summarise.inform = FALSE, width = 160)
set.seed(20261027)
NPERM <- 2000
hr    <- function(s) cat("\n", strrep("=", 78), "\n", s, "\n", strrep("=", 78), "\n", sep = "")
pr    <- function(d) print(as.data.frame(d), row.names = FALSE)
rnum  <- function(r) suppressWarnings(as.numeric(gsub("\\D", "", r)))
fmt_p <- function(p) ifelse(is.na(p), "p n/a",
           ifelse(p <= 1 / (NPERM + 1) + 1e-12, sprintf("p <= %.4f", 1 / (NPERM + 1)),
                  paste0("p = ", formatC(signif(p, 2), format = "g"))))

need <- c("DR/out/DR21_copynumber.csv", "DR/out/DR02_propagated_blocks.csv", "DR/out/pairwise_ks.csv")
miss <- need[!file.exists(need)]
if (length(miss)) stop("missing input: ", paste(miss, collapse = ", "),
                       " (run from the subgenome_prototype directory)")
dir.create("DR/fig", showWarnings = FALSE)

HEX    <- c("regia", "binata", "paradoxa", "scorpioides")
SPP    <- c(HEX, "capensis")
POOLED <- "4 hexaploids pooled"
CONST  <- list(regia = c(2, 1), binata = c(2, 1), paradoxa = c(2, 1),
               scorpioides = c(2, 1), capensis = c(4, 2))      # A sets, B sets
comp_levels <- function(k) vapply(0:k, function(a) paste0(strrep("A", a), strrep("B", k - a)), "")
hyper <- function(k, nA, nB) {   # composition if k of the copies survive at random
  a <- 0:k
  p <- if (k > nA + nB) rep(0, k + 1) else choose(nA, a) * choose(nB, k - a) / choose(nA + nB, k)
  setNames(p, comp_levels(k))
}
FILL3 <- c(BBB = "#993C1D", ABB = "#F0997B", AAB = "#1D9E75", AAA = "#9FE1CB")

# =============================================================================
hr("0. RECOVER THE FIGURE - DR21 loci with exactly three labelled copies")
CN <- read_csv("DR/out/DR21_copynumber.csv", show_col_types = FALSE)
stopifnot(all(c("genome", "locus", "A", "B", "sp") %in% names(CN)))
unk <- setdiff(unique(CN$sp), SPP)
if (length(unk)) stop("species not set up in this script: ", paste(unk, collapse = ", "))
CN <- CN %>% mutate(ncp = A + B,
                    region = ifelse(grepl("-", locus), sub("-.*$", "", locus), NA_character_))
reg <- unique(CN$region[!is.na(CN$region)]); lev <- reg[order(rnum(reg))]

binom_table <- function(d, kk, set) {
  d <- d %>% filter(ncp == kk)
  if (!nrow(d)) return(NULL)
  pA  <- sum(d$A) / (kk * nrow(d))
  obs <- tabulate(d$A + 1, nbins = kk + 1)
  ex  <- nrow(d) * dbinom(0:kk, kk, pA)
  tibble(set = set, k = kk, loci = nrow(d), pA = round(pA, 3), comp = comp_levels(kk),
         observed = obs, independent = round(ex, 1), ratio = round(obs / ex, 2))
}
R0 <- bind_rows(binom_table(CN, 3, "all five species"),
                binom_table(CN %>% filter(sp %in% HEX), 3, "four hexaploids"))
pr(R0)
cat("\n  The figure: 951 loci; BBB 0.2x, ABB 0.5x, AAB 1.7x, AAA 0.4x. If the 'all five\n")
cat("  species' rows give that, the figure is DR21's A + B = 3 loci pooled over all\n")
cat("  five species, capensis included.\n")
write_csv(R0, "DR/out/DR21c_0_recovered.csv")

# =============================================================================
hr("1. PER SPECIES - observed against independent labels at the species' own A-rate")
R1 <- bind_rows(lapply(SPP, function(s) bind_rows(lapply(2:3, function(kk)
  binom_table(CN %>% filter(sp == s), kk, s)))))
R1$constitution <- round(unname(mapply(function(s, kk, comp, loci) {
  cc <- CONST[[s]]; loci * hyper(kk, cc[1], cc[2])[[comp]] },
  R1$set, R1$k, R1$comp, R1$loci)), 1)
cat("  constitution = expected count if each locus kept a random subset of AAB's\n")
cat("  copies (capensis AAAABB). For a hexaploid every 3-copy locus is AAB and every\n")
cat("  2-copy locus is AA or AB.\n\n  3-copy loci\n")
pr(R1 %>% filter(k == 3))
cat("\n  2-copy loci\n")
pr(R1 %>% filter(k == 2))
V1 <- R1 %>% filter(k == 3, set %in% HEX) %>% select(set, loci, comp, observed) %>%
  pivot_wider(names_from = comp, values_from = observed, values_fill = 0) %>%
  mutate(not_AAB_pct = round(100 * (AAA + ABB + BBB) / loci, 1),
         B_called_A  = round(100 * AAA / loci, 1),
         A_called_B  = round(100 * ABB / loci / 2, 1))
cat("\n  hexaploid 3-copy loci that break AAB. One wrong copy turns AAB into AAA (a B\n")
cat("  copy called A, or a tandem A) or ABB (one of the two A copies called B, or a\n")
cat("  tandem B): per-copy rates ~ AAA share and ABB share / 2. Upper bounds.\n\n")
pr(V1)
c6 <- CN %>% filter(sp == "capensis", ncp == 6) %>% count(A, B, name = "loci")
cat("\n  capensis loci with all six copies (AAAABB predicts A = 4, B = 2):\n")
if (nrow(c6)) pr(c6) else cat("  none\n")
write_csv(R1, "DR/out/DR21c_1_per_species.csv")

# =============================================================================
hr("2. HONEST NULL - shuffle segment labels within each species")
PB <- read_csv("DR/out/DR02_propagated_blocks.csv", show_col_types = FALSE) %>%
  distinct(genome, gene, .keep_all = TRUE)
stopifnot(all(c("genome", "gene", "chr", "region", "label") %in% names(PB)))
KS <- read_csv("DR/out/pairwise_ks.csv", col_types = cols(.default = col_character()))
stopifnot(all(c("anchor", "seq1", "seq2", "sp1", "sp2") %in% names(KS)))
DROS <- unique(CN$genome)
TIPS <- bind_rows(KS %>% transmute(anchor, genome = sp1, seq = seq1),
                  KS %>% transmute(anchor, genome = sp2, seq = seq2)) %>%
  filter(genome %in% DROS) %>% distinct() %>% mutate(gene = sub("^[^@]*@", "", seq))
rm(KS)
# gene IDs may differ between files by isoform suffix or a 'both_agree_' prefix:
# per genome, the strictest key that matches best without creating duplicates
KEYS <- list(exact      = function(g) g,
             no_isoform = function(g) sub("\\.t[0-9]+$", "", g),
             no_prefix  = function(g) sub("^both_agree_", "", sub("\\.t[0-9]+$", "", g)))
pick <- bind_rows(lapply(DROS, function(g) {
  pg <- PB$gene[PB$genome == g]; tg <- TIPS$gene[TIPS$genome == g]
  rates <- sapply(KEYS, function(f) { kp <- f(pg)
    if (anyDuplicated(kp)) NA_real_ else mean(f(tg) %in% kp) })
  data.frame(genome = g, key = if (all(is.na(rates))) "exact" else names(which.max(rates)),
             stringsAsFactors = FALSE)
}))
key_of <- function(genome, gene) {
  out <- gene
  for (g in intersect(unique(genome), pick$genome)) {
    i <- genome == g; out[i] <- KEYS[[pick$key[pick$genome == g]]](gene[i])
  }
  out
}
PB   <- PB   %>% mutate(key = key_of(genome, gene))
TIPS <- TIPS %>% mutate(key = key_of(genome, gene))
TL <- TIPS %>%
  left_join(PB %>% select(genome, key, chr, region_pb = region, label), by = c("genome", "key")) %>%
  mutate(lab = ifelse(label %in% c("A", "B"), label, NA_character_))
REC <- TL %>% filter(!is.na(lab)) %>% count(genome, anchor, lab) %>%
  pivot_wider(names_from = lab, values_from = n, values_fill = 0)
for (s in c("A", "B")) if (!s %in% names(REC)) REC[[s]] <- 0L
CMP <- full_join(CN %>% select(genome, locus, A, B),
                 REC %>% transmute(genome, locus = anchor, A_r = A, B_r = B),
                 by = c("genome", "locus")) %>%
  mutate(status = case_when(is.na(A_r) ~ "only_DR21", is.na(A) ~ "only_rebuilt",
                            A == A_r & B == B_r ~ "identical", TRUE ~ "different"))
OK <- CMP %>% filter(status == "identical") %>% select(genome, locus)
repro <- CN %>% filter(ncp %in% 2:3) %>%
  left_join(CMP %>% select(genome, locus, status), by = c("genome", "locus")) %>%
  group_by(sp, copies = ncp) %>%
  summarise(loci = n(), rebuilt_exactly = sum(status == "identical", na.rm = TRUE)) %>%
  ungroup() %>% mutate(pct = round(100 * rebuilt_exactly / loci, 1)) %>%
  left_join(pick %>% transmute(sp = sub("^Drosera_", "", genome), key), by = "sp")
cat("  DR21 loci whose tips were rebuilt exactly (pairwise_ks tips + DR02 labels);\n")
cat("  the shuffle runs on these. key = how gene IDs were matched between files.\n\n")
pr(repro)

perm_one <- function(s, kk) {
  g   <- unique(CN$genome[CN$sp == s])
  ids <- CN %>% filter(genome == g, ncp == kk) %>% semi_join(OK, by = c("genome", "locus")) %>%
    pull(locus)
  if (length(ids) < 10) return(NULL)
  U  <- PB %>% filter(genome == g, label %in% c("A", "B")) %>% distinct(region, chr, label) %>%
    mutate(unit = paste(region, chr, label, sep = "|"))
  tp <- TL %>% filter(genome == g, anchor %in% ids, !is.na(lab)) %>%
    mutate(unit = paste(region_pb, chr, lab, sep = "|")) %>% arrange(anchor, unit)
  if (!all(table(tp$anchor) == kk)) return(NULL)
  idx <- match(tp$unit, U$unit)
  if (anyNA(idx)) return(NULL)
  M   <- matrix(idx, ncol = kk, byrow = TRUE)
  cnt <- function(lbl) tabulate(rowSums(matrix(lbl[M] == "A", ncol = kk)) + 1, nbins = kk + 1)
  obs <- cnt(U$label)
  dr  <- CN %>% filter(genome == g, locus %in% ids)
  if (!identical(as.integer(obs), tabulate(dr$A + 1, nbins = kk + 1)))
    cat("  WARNING: rebuilt compositions differ from DR21 for", s, kk, "copies\n")
  list(obs = obs, PM = t(replicate(NPERM, cnt(sample(U$label)))), loci = length(ids),
       segments = nrow(U), segments_A = sum(U$label == "A"),
       tandem = mean(apply(M, 1, function(r) anyDuplicated(r) > 0)))
}
PERM <- list()
for (s in SPP) for (kk in 2:3) PERM[[paste(s, kk)]] <- perm_one(s, kk)

summ <- function(panel, kk, obs, PM, cons) {
  lev <- comp_levels(kk); n <- sum(obs)
  pA <- sum(obs * (0:kk)) / (kk * n)
  binom <- n * dbinom(0:kk, kk, pA)
  if (is.null(PM)) {
    lo <- hi <- p_hi <- p_lo <- rep(NA_real_, kk + 1); pD <- NA_real_
  } else {
    mu   <- colMeans(PM)
    lo   <- apply(PM, 2, quantile, 0.025, names = FALSE)
    hi   <- apply(PM, 2, quantile, 0.975, names = FALSE)
    p_hi <- (colSums(sweep(PM, 2, obs, ">=")) + 1) / (nrow(PM) + 1)
    p_lo <- (colSums(sweep(PM, 2, obs, "<=")) + 1) / (nrow(PM) + 1)
    Dfun <- function(x) sum(ifelse(mu > 0, (x - mu)^2 / mu, 0))
    pD   <- (sum(apply(PM, 1, Dfun) >= Dfun(obs)) + 1) / (nrow(PM) + 1)
  }
  tibble(panel = panel, k = kk, loci = n, pA = pA, comp = factor(lev, levels = lev),
         observed = as.numeric(obs), binom = binom, lo = lo, hi = hi,
         p_hi = p_hi, p_lo = p_lo, p_D = pD, constitution = n * unname(cons))
}
S <- list()
for (s in SPP) for (kk in 2:3) {
  cc <- CONST[[s]]; x <- PERM[[paste(s, kk)]]
  if (is.null(x)) {
    d <- CN %>% filter(sp == s, ncp == kk)
    if (nrow(d)) S[[paste(s, kk)]] <- summ(s, kk, tabulate(d$A + 1, nbins = kk + 1), NULL,
                                           hyper(kk, cc[1], cc[2]))
  } else S[[paste(s, kk)]] <- summ(s, kk, x$obs, x$PM, hyper(kk, cc[1], cc[2]))
}
for (kk in 2:3) {
  xs <- lapply(HEX, function(s) PERM[[paste(s, kk)]])
  if (!any(vapply(xs, is.null, TRUE))) {
    S[[paste(POOLED, kk)]] <- summ(POOLED, kk, Reduce(`+`, lapply(xs, `[[`, "obs")),
                                   Reduce(`+`, lapply(xs, `[[`, "PM")), hyper(kk, 2, 1))
  } else {
    d <- CN %>% filter(sp %in% HEX, ncp == kk)
    S[[paste(POOLED, kk)]] <- summ(POOLED, kk, tabulate(d$A + 1, nbins = kk + 1), NULL,
                                   hyper(kk, 2, 1))
  }
}
S <- bind_rows(S) %>% mutate(panel = factor(panel, levels = c(SPP, POOLED)))
cat("\n  observed against the shuffle (lo-hi = 95% of shuffles; p_excess / p_deficit =\n")
cat("  share of shuffles at least as high / as low as observed)\n\n")
pr(S %>% transmute(panel, k, loci, comp, observed, independent = round(binom, 1),
                   shuffle_lo = lo, shuffle_hi = hi,
                   p_excess = signif(p_hi, 2), p_deficit = signif(p_lo, 2)))
INFO <- bind_rows(lapply(names(PERM), function(nm) { x <- PERM[[nm]]
  tibble(species_copies = nm, loci = x$loci, segments = x$segments, segments_A = x$segments_A,
         two_on_one_segment_pct = round(100 * x$tandem, 1)) }))
cat("\n  segments shuffled, and loci with two copies on one segment (tandem-like):\n\n")
if (nrow(INFO)) pr(INFO) else cat("  shuffle not run: rebuilt tips did not reproduce DR21\n")
write_csv(S, "DR/out/DR21c_2_shuffle.csv")

# =============================================================================
hr("3. SENSITIVITY - drop 3-copy loci that also carry an unlabelled tip")
UNL <- TL %>% group_by(genome, anchor) %>% summarise(unlab = sum(is.na(lab))) %>% ungroup()
SE <- CN %>% filter(ncp == 3) %>% left_join(UNL, by = c("genome", "locus" = "anchor")) %>%
  mutate(clean = coalesce(unlab, 0L) == 0) %>%
  group_by(sp) %>%
  summarise(loci = n(), AAB_pct = round(100 * mean(A == 2), 1), clean_loci = sum(clean),
            AAB_pct_clean = round(100 * mean(A[clean] == 2), 1)) %>% ungroup()
cat("  such a locus really holds 4+ copies, so its composition says less.\n\n")
pr(SE)

# =============================================================================
hr("4. BY REGION - where do the 3-copy loci that are not AAB sit?")
G3 <- CN %>% filter(ncp == 3, !is.na(region)) %>%
  mutate(comp = factor(comp_levels(3)[A + 1], levels = comp_levels(3))) %>%
  count(sp, region, comp, name = "loci") %>%
  group_by(sp, region) %>% mutate(N = sum(loci), share = loci / N) %>% ungroup()
W4 <- G3 %>% group_by(sp, region) %>%
  summarise(cell = sprintf("%d/%d%s", sum(loci[comp == "AAB"]), first(N),
                           ifelse(sum(loci[comp == "ABB"]) > sum(loci[comp == "AAB"]), " <ABB", ""))) %>%
  ungroup() %>% pivot_wider(names_from = sp, values_from = cell) %>%
  select(region, any_of(SPP)) %>% mutate(region = factor(region, levels = lev)) %>% arrange(region)
cat("  cell = AAB loci / all 3-copy loci in that region; '<ABB' = more ABB than AAB.\n")
cat("  Marked in every hexaploid = suspect that Dionaea pair's orientation.\n")
cat("  capensis: about 60% AAB is what AAAABB with random loss gives.\n\n")
pr(W4)
write_csv(G3, "DR/out/DR21c_4_by_region.csv")

# =============================================================================
hr("5. FIGURES")
draw_k <- function(kk, file, title, subtitle) {
  d <- S %>% filter(k == kk)
  if (!nrow(d)) return(NULL)
  key <- if (kk == 3) "AAB" else "BB"
  tt <- d %>% group_by(panel) %>%
    summarise(loci = first(loci), obs_key = observed[comp == key], ch_key = binom[comp == key],
              p_key = if (kk == 3) p_hi[comp == key] else p_lo[comp == key], pD = first(p_D)) %>%
    ungroup() %>%
    mutate(lbl = ifelse(as.character(panel) == "capensis",
                        sprintf("capensis   n = %d\nnot independent: %s", loci, fmt_p(pD)),
                        sprintf("%s   n = %d\n%s %.1f%% vs %.1f%% by chance, %s", as.character(panel),
                                loci, key, 100 * obs_key / loci, 100 * ch_key / loci, fmt_p(p_key))))
  long <- d %>% select(panel, comp, observed, binom) %>%
    pivot_longer(c(observed, binom), names_to = "what", values_to = "count") %>%
    mutate(what = factor(ifelse(what == "observed", "observed", "if labels were independent"),
                         levels = c("observed", "if labels were independent")))
  g <- ggplot(long, aes(comp, count, fill = what)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.78) +
    geom_errorbar(data = d, aes(x = comp, ymin = lo, ymax = hi), inherit.aes = FALSE,
                  position = position_nudge(x = 0.2), width = 0.18, colour = "#444441") +
    geom_text(data = d, aes(x = comp, y = observed, label = sprintf("%.1fx", observed / binom)),
              inherit.aes = FALSE, position = position_nudge(x = -0.2), vjust = -0.5,
              size = 2.9, colour = "#0F6E56", fontface = "bold") +
    facet_wrap(~ panel, nrow = 2, scales = "free_y",
               labeller = as_labeller(setNames(tt$lbl, as.character(tt$panel)))) +
    scale_fill_manual(values = c("observed" = "#1D9E75", "if labels were independent" = "#B4B2A9"),
                      name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(title = title, subtitle = subtitle, x = "labels of the copies at one locus", y = "loci") +
    theme_minimal(10) + theme(legend.position = "top", strip.text = element_text(hjust = 0))
  if (kk == 3) g <- g + geom_point(data = d, aes(x = comp, y = constitution), inherit.aes = FALSE,
                                   shape = 23, size = 2.3, fill = "black",
                                   position = position_nudge(x = -0.2))
  ggsave(file, g, width = 12, height = 7.5, dpi = 150, bg = "white")
  g
}
g1 <- draw_k(3, "DR/fig/DR21c_1_three_copy.png",
             "DR21c-1  Three-copy loci: the labels at one locus are not independent",
             paste0("grey = independent labels at the species' A-rate; whisker = 95% of ", NPERM,
                    " segment-label shuffles; diamond = AAB (capensis: AAAABB with random loss)"))
g2 <- draw_k(2, "DR/fig/DR21c_2_two_copy.png",
             "DR21c-2  Two-copy loci: AAB allows AA and AB, never BB",
             paste0("grey = independent labels at the species' A-rate; whisker = 95% of ", NPERM,
                    " segment-label shuffles"))

G3p <- G3 %>% mutate(region = factor(region, levels = lev), sp = factor(sp, levels = SPP))
g3 <- ggplot(G3p, aes(region, share, fill = comp)) +
  geom_col(width = 0.75) +
  geom_text(data = distinct(G3p, sp, region, N), aes(region, 1.07, label = N),
            inherit.aes = FALSE, size = 2.6) +
  facet_wrap(~ sp, ncol = 1) +
  scale_fill_manual(values = FILL3, name = NULL, drop = FALSE) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0%", "50%", "100%")) +
  labs(title = "DR21c-3  Three-copy loci, region by region",
       subtitle = "a mis-oriented Dionaea pair would turn one region orange in every species; number = loci",
       x = "region (Nepenthes chromosome = one Dionaea homeolog pair)", y = "share of 3-copy loci") +
  theme_minimal(10) + theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("DR/fig/DR21c_3_by_region.png", g3, width = 8, height = 10, dpi = 150, bg = "white")

RAT <- R1 %>% filter(k == 3) %>%
  transmute(sp = set, A = as.numeric(nchar(gsub("B", "", comp))), B = 3 - A, ratio)
G4 <- CN %>% count(sp, A, B, name = "loci") %>%
  mutate(A = as.numeric(A), B = as.numeric(B)) %>%
  left_join(RAT, by = c("sp", "A", "B")) %>%
  filter(A <= 6, B <= 4) %>%
  mutate(sp = factor(sp, levels = SPP), on_diag = A + B == 3)
g4 <- ggplot(G4, aes(B, A)) +
  geom_abline(slope = -1, intercept = 3, linetype = "dashed", colour = "grey60") +
  geom_point(aes(size = loci), colour = "#888780", alpha = 0.5) +
  geom_point(data = filter(G4, on_diag), aes(size = loci), shape = 21, fill = NA,
             colour = "black", stroke = 0.9) +
  geom_text(data = filter(G4, on_diag), aes(label = sprintf("%.1fx", ratio)), nudge_x = 0.45,
            size = 3, colour = "#0F6E56", fontface = "bold") +
  facet_wrap(~ sp, nrow = 1) +
  scale_size_area(max_size = 13, name = "loci") +
  scale_x_continuous(breaks = 0:4) + scale_y_continuous(breaks = 0:6) +
  coord_cartesian(xlim = c(-0.5, 4.6), ylim = c(-0.5, 6.5)) +
  labs(title = "DR21c-4  The three-copy figure is the A + B = 3 diagonal of the bubble plot",
       subtitle = "ringed = loci that kept exactly three labelled copies; green = observed / independent-label expectation",
       x = "B copies at a locus", y = "A copies at a locus") +
  theme_minimal(10) + theme(legend.position = "right")
ggsave("DR/fig/DR21c_4_bubble_link.png", g4, width = 14, height = 4.6, dpi = 150, bg = "white")
pdf("DR/fig/DR21c_all.pdf", width = 12, height = 8)
for (g in Filter(Negate(is.null), list(g1, g2, g3, g4))) print(g)
invisible(dev.off())
cat("  wrote DR/fig/DR21c_{1_three_copy,2_two_copy,3_by_region,4_bubble_link}.png and DR21c_all.pdf\n")

# =============================================================================
hr("SUMMARY")
line_of <- function(d, kk) {
  if (!nrow(d)) return("n/a")
  if (as.character(d$panel[1]) == "capensis")
    return(sprintf("n = %d, labels not independent %s", d$loci[1], fmt_p(d$p_D[1])))
  key <- if (kk == 3) "AAB" else "BB"; r <- d[d$comp == key, ]
  sprintf("n = %d, %s %.1f%% vs %.1f%% by chance, %s", r$loci, key,
          100 * r$observed / r$loci, 100 * r$binom / r$loci, fmt_p(if (kk == 3) r$p_hi else r$p_lo))
}
for (s in c(SPP, POOLED)) {
  cat(sprintf("  %-20s 3-copy  %s\n", s, line_of(S %>% filter(panel == s, k == 3), 3)))
  cat(sprintf("  %-20s 2-copy  %s\n", "", line_of(S %>% filter(panel == s, k == 2), 2)))
}
cat("\n  hexaploid 3-copy loci breaking AAB (implied per-copy flips, upper bounds):\n")
for (i in seq_len(nrow(V1)))
  cat(sprintf("  %-12s %4.1f%% not AAB   B called A ~%4.1f%%   A called B ~%4.1f%%\n",
              V1$set[i], V1$not_AAB_pct[i], V1$B_called_A[i], V1$A_called_B[i]))
FLAG <- G3 %>% filter(sp %in% HEX) %>% group_by(sp, region) %>%
  summarise(abb_wins = sum(loci[comp == "ABB"]) > sum(loci[comp == "AAB"])) %>% ungroup() %>%
  group_by(region) %>% summarise(n_sp = n(), n_abb = sum(abb_wins)) %>% ungroup()
fl_all <- FLAG$region[FLAG$n_abb == FLAG$n_sp]
fl_any <- FLAG$region[FLAG$n_abb > 0 & FLAG$n_abb < FLAG$n_sp]
cat(sprintf("\n  regions with more ABB than AAB in every hexaploid: %s\n",
            if (length(fl_all)) paste(fl_all, collapse = ", ") else "none"))
cat(sprintf("  ... in some hexaploids only: %s\n",
            if (length(fl_any)) paste(fl_any, collapse = ", ") else "none"))
