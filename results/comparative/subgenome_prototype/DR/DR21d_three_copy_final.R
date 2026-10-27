#!/usr/bin/env Rscript
# =============================================================================
# DR21d - FINAL THREE-COPY FIGURE (the original design), REBUILT FROM DR21
# =============================================================================
# At loci where a species kept exactly three labelled copies: how the three are
# labelled, against independent labels at the pooled A-rate. Same design as the
# original figure, but every number comes from DR/out/DR21_copynumber.csv
# (copies per locus by label), and only the four hexaploids enter: capensis is
# AAAABB, so its 3-copy loci are 3 of its 6 copies, with another expectation.
#
# The p in the subtitle is DR21c's segment-label shuffle (2000 shuffles). Labels
# are set per segment, so loci in one region are not independent; the shuffle
# keeps that, the binomial does not. The binomial bars stay as the picture.
#
# Also prints the caption numbers: per species, and for paradoxa where its
# non-AAB loci sit and how old their same-label pairs are (one homeolog recently
# replaced by a copy of another, or a mislabelled segment?).
#
# Inputs  DR/out/DR21_copynumber.csv; DR21c_2_shuffle.csv if present; for the
#         paradoxa pairs DR02_propagated_blocks.csv and pairwise_ks.csv
# Outputs DR/fig/DR21d_three_copy.{png,pdf}; DR/out/DR21d_*.csv
# =============================================================================

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr) })
options(dplyr.summarise.inform = FALSE, width = 150)
hr <- function(s) cat("\n", strrep("=", 78), "\n", s, "\n", strrep("=", 78), "\n", sep = "")
pr <- function(d) print(as.data.frame(d), row.names = FALSE)
if (!file.exists("DR/out/DR21_copynumber.csv"))
  stop("DR/out/DR21_copynumber.csv not found (run from the subgenome_prototype directory)")
dir.create("DR/fig", showWarnings = FALSE)

HEX   <- c("regia", "binata", "paradoxa", "scorpioides")
COMP  <- c("BBB", "ABB", "AAB", "AAA")          # 0, 1, 2, 3 copies labelled A
GREEN <- "#1D9E75"; GREY <- "#B4B2A9"; DARK <- "#0F6E56"

CN <- read_csv("DR/out/DR21_copynumber.csv", show_col_types = FALSE) %>%
  mutate(ncp = A + B, region = ifelse(grepl("-", locus), sub("-.*$", "", locus), NA_character_))
T3  <- CN %>% filter(sp %in% HEX, ncp == 3)
NL  <- nrow(T3)
PA  <- sum(T3$A) / (3 * NL)
OBS <- setNames(tabulate(T3$A + 1, nbins = 4), COMP)
IND <- setNames(NL * dbinom(0:3, 3, PA), COMP)

# honest p: DR21c's segment-label shuffle for the pooled hexaploids
SH <- NULL; p_txt <- NA_character_
if (file.exists("DR/out/DR21c_2_shuffle.csv")) {
  SH  <- read_csv("DR/out/DR21c_2_shuffle.csv", show_col_types = FALSE)
  row <- SH %>% filter(panel == "4 hexaploids pooled", k == 3, comp == "AAB")
  if (nrow(row) == 1 && !is.na(row$p_hi)) {
    p_txt <- if (row$p_hi < 0.001) "p < 0.001" else sprintf("p = %.2g", row$p_hi)
    if (row$loci != NL)
      cat(sprintf("  NOTE: DR21c's shuffle ran on %d of these %d loci (those it rebuilt exactly)\n",
                  as.integer(row$loci), NL))
  }
}
SUB <- paste0(format(NL, big.mark = ","), " three-copy loci in four hexaploid species",
              if (is.na(p_txt)) "" else paste0("  |  shuffled segment labels: ", p_txt))

draw <- function() {
  op <- par(mar = c(5.2, 6, 7, 2), xpd = NA)
  M  <- rbind(OBS, IND)
  bp <- barplot(M, beside = TRUE, space = c(0, 0.7), col = c(GREEN, GREY), border = NA,
                ylim = c(0, max(M) * 1.15), axes = FALSE, names.arg = COMP, cex.names = 1.35)
  axis(2, las = 1, cex.axis = 1.1)
  text(bp[1, ], OBS, sprintf("%.1fx", OBS / IND), pos = 3, col = DARK, font = 2, cex = 1.15)
  mtext("three-copy loci", side = 2, line = 4.2, cex = 1.25)
  mtext("how the three copies at one locus are labelled", side = 1, line = 3.6, cex = 1.25)
  mtext("Copies at the same locus are not labelled independently", side = 3, line = 4.6,
        adj = 0, font = 2, cex = 1.55)
  mtext(SUB, side = 3, line = 3.0, adj = 0, col = "grey40", cex = 1.12)
  legend("topleft", legend = c("observed", "if the labels were independent"),
         fill = c(GREEN, GREY), border = NA, bty = "n", cex = 1.2)
  par(op)
}
png("DR/fig/DR21d_three_copy.png", width = 1420, height = 880, res = 150); draw(); invisible(dev.off())
pdf("DR/fig/DR21d_three_copy.pdf", width = 9.5, height = 5.9); draw(); invisible(dev.off())

# =============================================================================
hr("NUMBERS FOR THE FIGURE AND THE CAPTION")
FIG <- tibble(labels = COMP, observed = as.integer(OBS), independent = round(IND, 1),
              ratio = round(OBS / IND, 2))
cat(sprintf("  four hexaploids pooled: %d loci, A-rate %.3f; AAB %.1f%% observed vs %.1f%% if independent\n\n",
            NL, PA, 100 * OBS[["AAB"]] / NL, 100 * IND[["AAB"]] / NL))
pr(FIG)
PS <- T3 %>% group_by(sp) %>%
  summarise(loci = n(), pA = sum(A) / (3 * n()), AAB = sum(A == 2)) %>% ungroup() %>%
  mutate(AAB_pct = round(100 * AAB / loci, 1), chance_pct = round(100 * 3 * pA^2 * (1 - pA), 1)) %>%
  select(sp, loci, AAB_pct, chance_pct)
if (!is.null(SH)) PS <- PS %>%
  left_join(SH %>% filter(k == 3, comp == "AAB") %>%
              transmute(sp = as.character(panel), shuffle_p = signif(p_hi, 2)), by = "sp")
cat("\n  per species (shuffle p from DR21c)\n\n")
pr(PS %>% mutate(sp = factor(sp, levels = HEX)) %>% arrange(sp))
write_csv(FIG, "DR/out/DR21d_figure_numbers.csv")
write_csv(PS,  "DR/out/DR21d_per_species.csv")

# =============================================================================
hr("PARADOXA - where its non-AAB loci sit, and how old their same-label pairs are")
P3 <- T3 %>% filter(sp == "paradoxa", !is.na(region)) %>% mutate(lab3 = COMP[A + 1]) %>%
  count(region, lab3, name = "loci") %>%
  pivot_wider(names_from = lab3, values_from = loci, values_fill = 0)
for (cc in COMP) if (!cc %in% names(P3)) P3[[cc]] <- 0L
P3 <- P3 %>% mutate(loci = BBB + ABB + AAB + AAA, AAB_pct = round(100 * AAB / loci, 1)) %>%
  select(region, loci, all_of(COMP), AAB_pct) %>% arrange(desc(loci))
pr(P3)
flag <- P3$region[P3$AAB < P3$loci - P3$AAB]
oth  <- P3 %>% filter(!region %in% flag)
cat(sprintf("\n  regions where non-AAB loci outnumber AAB: %s\n",
            if (length(flag)) paste(flag, collapse = ", ") else "none"))
cat(sprintf("  all other regions together: %d of %d loci AAB (%.1f%%)\n",
            as.integer(sum(oth$AAB)), as.integer(sum(oth$loci)),
            100 * sum(oth$AAB) / max(1, sum(oth$loci))))

need2 <- c("DR/out/DR02_propagated_blocks.csv", "DR/out/pairwise_ks.csv")
if (all(file.exists(need2))) {
  PB <- read_csv(need2[1], show_col_types = FALSE) %>% filter(genome == "Drosera_paradoxa") %>%
    distinct(gene, .keep_all = TRUE)
  KS <- read_csv(need2[2], col_types = cols(.default = col_character())) %>%
    filter(sp1 == "Drosera_paradoxa", sp2 == "Drosera_paradoxa") %>%
    mutate(dS = suppressWarnings(as.numeric(dS)),
           g1 = sub("^[^@]*@", "", seq1), g2 = sub("^[^@]*@", "", seq2))
  # gene IDs matched as in DR21c: the strictest key that matches best without duplicates
  KEYS <- list(exact      = function(g) g,
               no_isoform = function(g) sub("\\.t[0-9]+$", "", g),
               no_prefix  = function(g) sub("^both_agree_", "", sub("\\.t[0-9]+$", "", g)))
  rate  <- sapply(KEYS, function(f) { kp <- f(PB$gene)
    if (anyDuplicated(kp)) NA_real_ else mean(f(c(KS$g1, KS$g2)) %in% kp) })
  kname <- if (all(is.na(rate))) "exact" else names(which.max(rate))
  kf    <- KEYS[[kname]]
  LAB   <- setNames(PB$label, kf(PB$gene))
  L3 <- T3 %>% filter(sp == "paradoxa") %>% transmute(anchor = locus, lab3 = COMP[A + 1], region)
  PP <- KS %>% mutate(l1 = unname(LAB[kf(g1)]), l2 = unname(LAB[kf(g2)])) %>%
    filter(l1 %in% c("A", "B"), l2 %in% c("A", "B"), !is.na(dS), dS >= 0, dS < 5) %>%
    inner_join(L3, by = "anchor") %>%
    mutate(pair = ifelse(l1 == l2, paste0(l1, l2), "AB"))
  cat(sprintf("\n  gene IDs matched with key '%s' (%.0f%% of pair members found in DR02)\n", kname,
              100 * ifelse(is.na(rate[[kname]]), 0, rate[[kname]])))
  ref <- PP %>% filter(lab3 == "AAB") %>% group_by(pair) %>%
    summarise(pairs = n(), median_dS = round(median(dS), 3)) %>% ungroup()
  ab  <- ref$median_dS[ref$pair == "AB"]; ab <- if (length(ab)) ab else NA_real_
  cat("  reference depths at paradoxa's AAB loci: the AA pair is A vs A', the AB pairs A vs B\n\n")
  pr(ref)
  YG <- PP %>% filter(lab3 != "AAB", pair %in% c("AA", "BB")) %>%
    group_by(anchor, region, lab3) %>% summarise(youngest = min(dS)) %>% ungroup() %>%
    group_by(region, labels = lab3) %>%
    summarise(loci = n(), median_youngest_dS = round(median(youngest), 3),
              young_pct = round(100 * mean(youngest < 0.3 * ab), 1)) %>% ungroup() %>%
    arrange(desc(loci))
  cat("\n  non-AAB loci: dS of the youngest same-label pair; young = below 0.3 x the A vs B depth\n")
  cat("  young      -> one copy recently replaced by a copy of its homeolog (conversion,\n")
  cat("                exchange) or a recent duplicate: lineage-specific, AAB untouched\n")
  cat("  ~ A vs A'  -> the same-label pair is really A and A': a mislabelled segment\n\n")
  pr(YG)
  write_csv(YG, "DR/out/DR21d_paradoxa_pairs.csv")
} else cat("\n  pair ages skipped: needs DR02_propagated_blocks.csv and pairwise_ks.csv\n")

cat("\n  figure: DR/fig/DR21d_three_copy.png (+ .pdf)\n")
