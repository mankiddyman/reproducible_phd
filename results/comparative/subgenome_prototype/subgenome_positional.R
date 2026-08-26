#!/usr/bin/env Rscript
# PIN LIBPATH: stale ~/R user library shadows the env lib (rlang 1.1.3 vs dplyr needing >=1.1.7)
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Where does the dS asymmetry live: tract, chromosome-wide, or nowhere?
# sg1/sg2 is arbitrary BETWEEN pairs by construction. Consistent WITHIN a pair,
# which is all the paired test needs. Global orientation is the unknown.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ragg)
})

WD  <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype"
GS  <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results/combBed.txt"
FAI <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/Dionaea_muscipula/assembly_final/external_collapsed/Dionaea_muscipula_chr.fa.fai"
setwd(WD)

DS_MAX <- 1.5     # censor saturated / misaligned
DS_CONV <- 0.10   # below this, homeologs likely gene-converted
WIN_K <- 25; WIN_S <- 5

d <- read_csv("distances.csv", show_col_types = FALSE)

bed <- read.table(GS, header = TRUE, sep = "\t", quote = "", comment.char = "",
                  stringsAsFactors = FALSE)
message("combBed columns: ", paste(names(bed), collapse = ", "))

pick <- function(nms, cands) { h <- cands[cands %in% nms]; if (!length(h)) NA_character_ else h[1] }
CHR <- pick(names(bed), c("chr", "chromosome"))
ST  <- pick(names(bed), c("start", "st", "begin"))
GEN <- pick(names(bed), c("genome", "genomeID", "gid"))
if (any(is.na(c(CHR, ST, GEN)))) stop("combBed missing chr/start/genome")

dio_ids <- unique(c(d$dio1, d$dio2))
cands <- intersect(c("id", "ofID", "gene", "geneID", "name"), names(bed))
hr <- sapply(cands, function(cc) mean(dio_ids %in% bed[[cc]]))
message("id-column match rates: ", paste(sprintf("%s=%.3f", cands, hr), collapse = "  "))
if (!length(cands) || max(hr) < 0.5) stop("no combBed id column matches distances.csv")
ID <- cands[which.max(hr)]

pos <- bed %>% filter(.data[[GEN]] == "Dionaea_muscipula") %>%
  transmute(gene = .data[[ID]], p_start = as.numeric(.data[[ST]])) %>%
  distinct(gene, .keep_all = TRUE)

## ---- keep everything; flag same-chr rather than dropping --------------------
d2 <- d %>%
  filter(is.finite(dS_nep_dio1), is.finite(dS_nep_dio2), is.finite(dS_dio1_dio2)) %>%
  mutate(same_chr   = chr1 == chr2,
         pl1        = sub("_sg[0-9]+_s[0-9]+$", "", chr1),
         pl2        = sub("_sg[0-9]+_s[0-9]+$", "", chr2),
         cross_pair = pl1 != pl2,
         pair_lab   = ifelse(cross_pair, NA_character_, pl1))

message(sprintf("triplets %d | same-chr (HE) %d | cross-pair (inter-pair translocation) %d | within-pair %d",
                nrow(d2), sum(d2$same_chr), sum(d2$cross_pair),
                sum(!d2$same_chr & !d2$cross_pair)))

## ---- paired test set: orientation consistent WITHIN pair only ---------------
ph <- d2 %>% filter(!same_chr, !cross_pair) %>%
  mutate(swap  = chr2 < chr1,
         chrA  = pmin(chr1, chr2), chrB = pmax(chr1, chr2),
         geneA = ifelse(swap, dio2, dio1),
         dS_A  = ifelse(swap, dS_nep_dio2, dS_nep_dio1),
         dS_B  = ifelse(swap, dS_nep_dio1, dS_nep_dio2),
         diff  = dS_A - dS_B)

raw  <- filter(ph, dS_A < DS_MAX, dS_B < DS_MAX)
conv <- filter(raw, dS_dio1_dio2 >= DS_CONV)
message(sprintf("censor dS>%.1f: %d -> %d | conversion filter dS_hom<%.2f: %d -> %d",
                DS_MAX, nrow(ph), nrow(raw), DS_CONV, nrow(raw), nrow(conv)))

conv <- left_join(conv, pos, by = c("geneA" = "gene"))
message(sprintf("position join: %.1f%% placed", 100 * mean(!is.na(conv$p_start))))

hl <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 6) return(tibble(n = length(x), est = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
  w <- suppressWarnings(wilcox.test(x, conf.int = TRUE))
  tibble(n = length(x), est = unname(w$estimate), lo = w$conf.int[1], hi = w$conf.int[2], p = w$p.value)
}

stats_for <- function(df, tag) {
  df %>% group_by(pair_lab) %>%
    reframe(chrA = sort(unique(chrA))[1], chrB = sort(unique(chrB))[1], hl(diff)) %>%
    mutate(set = tag, p_adj = p.adjust(p, "BH"),
           faster = ifelse(is.na(est), NA, ifelse(est > 0, chrA, chrB)))
}
st <- bind_rows(stats_for(raw, "censored_only"), stats_for(conv, "conv_filtered"))
write_csv(st, "phasing_by_chrpair_v2.csv")

cochran <- function(s) {
  s <- filter(s, is.finite(est), is.finite(lo), hi > lo)
  se <- (s$hi - s$lo) / (2 * qnorm(0.975)); w <- 1 / se^2
  mu <- sum(w * s$est) / sum(w); Q <- sum(w * (s$est - mu)^2); df <- nrow(s) - 1
  tibble(pooled = mu, pooled_se = sqrt(1/sum(w)), Q = Q, df = df,
         p_Q = pchisq(Q, df, lower.tail = FALSE), I2 = max(0, (Q - df)/Q))
}
het <- bind_rows(cochran(filter(st, set == "censored_only")) %>% mutate(set = "censored_only"),
                 cochran(filter(st, set == "conv_filtered")) %>% mutate(set = "conv_filtered"))
write_csv(het, "phasing_heterogeneity.csv")

runs_test <- function(s) {
  s <- s[s != 0 & !is.na(s)]; n1 <- sum(s > 0); n2 <- sum(s < 0); n <- n1 + n2
  if (n1 < 5 || n2 < 5) return(tibble(runs = NA_real_, exp_runs = NA_real_, z = NA_real_, p = NA_real_))
  r <- 1 + sum(s[-1] != s[-n]); mu <- 2*n1*n2/n + 1
  v <- 2*n1*n2*(2*n1*n2 - n)/(n^2*(n-1)); z <- (r - mu)/sqrt(v)
  tibble(runs = r, exp_runs = mu, z = z, p = 2*pnorm(-abs(z)))
}
ordered <- conv %>% filter(!is.na(p_start)) %>% arrange(pair_lab, p_start)
rt <- ordered %>% group_by(pair_lab) %>% reframe(runs_test(sign(diff))) %>%
  mutate(p_adj = p.adjust(p, "BH"))
write_csv(rt, "positional_runs.csv")

win <- ordered %>% group_by(pair_lab) %>% group_modify(function(g, key) {
  n <- nrow(g); if (n < WIN_K) return(tibble())
  s <- seq(1, n - WIN_K + 1, by = WIN_S)
  tibble(mid = sapply(s, function(i) median(g$p_start[i:(i+WIN_K-1)])),
         med = sapply(s, function(i) median(g$diff[i:(i+WIN_K-1)])))
}) %>% ungroup()

fai <- read.table(FAI, stringsAsFactors = FALSE)[,1:2]; names(fai) <- c("chr","len")
lenc <- st %>% filter(set == "conv_filtered") %>%
  left_join(fai, by = c("chrA"="chr")) %>% rename(lenA = len) %>%
  left_join(fai, by = c("chrB"="chr")) %>% rename(lenB = len) %>%
  mutate(dlen_Mb = (lenA - lenB)/1e6) %>%
  select(pair_lab, chrA, chrB, lenA, lenB, dlen_Mb, est, lo, hi, n)
write_csv(lenc, "phasing_length_covariate.csv")
ct <- with(filter(lenc, is.finite(est)), cor.test(dlen_Mb, est, method = "spearman"))

pE <- ggplot(ordered, aes(p_start/1e6, diff)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_point(alpha = 0.35, size = 0.9, colour = "#7570b3") +
  { if (nrow(win)) geom_line(data = win, aes(mid/1e6, med), colour = "#d95f02", linewidth = 0.9) } +
  facet_wrap(~ pair_lab, scales = "free_x", ncol = 2) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  labs(title = "E.  dS asymmetry along the chromosome",
       subtitle = paste0("Orange = ", WIN_K, "-gene sliding median. Flat at zero = no signal; uniform offset = chromosome-wide; localised block = exchange tract or assembly artifact."),
       x = "position on first chromosome of pair (Mb)", y = "dS(Nep, A) - dS(Nep, B)") +
  theme_bw(base_size = 11)

sf <- filter(st, set == "conv_filtered")
pF <- ggplot(sf, aes(est, reorder(pair_lab, est))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.22, orientation = "y", colour = "#1b9e77") +
  geom_point(aes(size = n), colour = "#1b9e77") +
  scale_size_continuous(range = c(1.5, 4)) +
  labs(title = "F.  Effect size per chromosome pair (Hodges-Lehmann, 95% CI)",
       subtitle = sprintf("pooled %.4f (SE %.4f) | Cochran Q p=%.3g, I2=%.0f%%",
                          het$pooled[het$set=="conv_filtered"], het$pooled_se[het$set=="conv_filtered"],
                          het$p_Q[het$set=="conv_filtered"], 100*het$I2[het$set=="conv_filtered"]),
       x = "median paired difference in dS", y = NULL, size = "gene pairs") +
  theme_bw(base_size = 11)

ggsave("subgenome_E_position.png", pE, width = 11, height = 9, dpi = 150, device = agg_png)
ggsave("subgenome_E_position.pdf", pE, width = 11, height = 9)
ggsave("subgenome_F_forest.png", pF, width = 8, height = 5, dpi = 150, device = agg_png)
ggsave("subgenome_F_forest.pdf", pF, width = 8, height = 5)

cat("\n=== effect sizes (conversion-filtered) ===\n"); print(as.data.frame(sf), digits = 3)
cat("\n=== heterogeneity ===\n"); print(as.data.frame(het), digits = 3)
cat("\n=== runs test (low z = same-sign clustering = tract) ===\n"); print(as.data.frame(rt), digits = 3)
cat("\n=== length covariate ===\n"); print(as.data.frame(lenc), digits = 3)
cat(sprintf("Spearman(dlen, effect) rho=%.3f p=%.3f\n", ct$estimate, ct$p.value))
