#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Per gene along a Dionaea chromosome pair: does each Drosera copy sit closer to the
# sgA or sgB Dionaea homeolog? Votes kept PER LINEAGE (species + chromosome of the
# voting gene, read directly from the tip, no join). No dominance used anywhere.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ape); library(ggplot2); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
WIN <- 15; MAXGAP <- 10e6; B <- 999
set.seed(1)

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
meta <- read_tsv("wgd7/tip_meta.tsv", show_col_types = FALSE) %>%
  left_join(bed %>% transmute(gene = id, genome, mid = (start + end)/2),
            by = c("gene", "genome"))
key <- function(x) gsub("@", "_", gsub("['\"]", "", x))

fs <- list.files("wgd7iq/tre", "\\.tre$", full.names = TRUE)
cat(sprintf("trees: %d\n", length(fs)))

one <- function(f) {
  tr <- tryCatch(read.tree(f), error = function(e) NULL)
  if (is.null(tr) || Ntip(tr) < 4) return(NULL)
  m <- meta[match(key(tr$tip.label), key(meta$tip)), ]
  if (any(is.na(m$genome))) return(NULL)
  nep <- tr$tip.label[m$genome == "Nepenthes_gracilis"]
  if (length(nep) != 1) return(NULL)
  tr <- tryCatch(root(tr, nep, resolve.root = TRUE), error = function(e) NULL)
  if (is.null(tr)) return(NULL)
  m <- meta[match(key(tr$tip.label), key(meta$tip)), ]

  di <- which(m$genome == "Dionaea_muscipula")
  if (length(di) != 2) return(NULL)
  c1 <- m$chr[di[1]]; c2 <- m$chr[di[2]]
  if (c1 == c2) return(NULL)
  p1 <- sub("_sg[0-9]+_s[0-9]+$", "", c1); p2 <- sub("_sg[0-9]+_s[0-9]+$", "", c2)
  if (p1 != p2) return(NULL)
  o <- order(c(c1, c2)); iA <- di[o[1]]; iB <- di[o[2]]
  tA <- tr$tip.label[iA]; tB <- tr$tip.label[iB]

  nd <- node.depth(tr)
  j <- which(!m$genome %in% c("Nepenthes_gracilis", "Dionaea_muscipula"))
  if (!length(j)) return(NULL)
  sup_all <- suppressWarnings(as.numeric(sub(".*/", "", tr$node.label)))
  med_sup <- if (all(is.na(sup_all))) NA_real_ else median(sup_all, na.rm = TRUE)
  do.call(rbind, lapply(j, function(k) {
    mA <- getMRCA(tr, c(tr$tip.label[k], tA)); mB <- getMRCA(tr, c(tr$tip.label[k], tB))
    dA <- nd[mA]; dB <- nd[mB]
    if (dA == dB) return(NULL)
    mv <- if (dA < dB) mA else mB
    sup <- suppressWarnings(as.numeric(sub(".*/", "", tr$node.label[mv - Ntip(tr)])))
    data.frame(pair = p1, chrA = m$chr[iA], chrB = m$chr[iB],
               posA = m$mid[iA], anchor = sub("\\.tre$", "", basename(f)),
               sp = m$genome[k],
               lin_chr = m$chr[k],                      # <- FIX: read from the tip
               sp_gene = m$gene[k],
               vote = ifelse(dA < dB, "A", "B"),
               support = sup, med_sup = med_sup, stringsAsFactors = FALSE)
  }))
}

v <- bind_rows(lapply(fs, one)) %>% filter(!is.na(posA), !is.na(lin_chr)) %>%
  mutate(lineage = paste(sub("Drosera_", "", sp), sub("_hap1$|_collapsed$", "", lin_chr)))
write_csv(v, "tract_votes7iq.csv")
cat("\n=== UFBoot support on the node that decides each vote ===\n")
cat(sprintf("votes with a numeric support value: %.1f%%\n", 100*mean(!is.na(v$support))))
print(summary(v$support))
cat(sprintf("support >=70: %d (%.0f%%) | >=95: %d (%.0f%%)\n",
            sum(v$support >= 70, na.rm=TRUE), 100*mean(v$support >= 70, na.rm=TRUE),
            sum(v$support >= 95, na.rm=TRUE), 100*mean(v$support >= 95, na.rm=TRUE)))
cat(sprintf("votes: %d | genes: %d | pairs: %d | lineages: %d\n",
            nrow(v), n_distinct(v$anchor), n_distinct(v$pair), n_distinct(v$lineage)))

## ---- gene coverage: is the chromosome actually scanned, or just its ends? ----
cov <- v %>% distinct(pair, anchor, posA) %>% arrange(pair, posA) %>%
  group_by(pair) %>%
  summarise(n_genes = n(), span_Mb = (max(posA)-min(posA))/1e6,
            median_gap_Mb = median(diff(posA))/1e6,
            max_gap_Mb = max(diff(posA))/1e6,
            frac_span_in_maxgap = max(diff(posA))/(max(posA)-min(posA)), .groups = "drop")
cat("\n=== gene coverage per chromosome pair ===\n")
print(as.data.frame(cov), digits = 3)
cat("large max_gap / high frac_span_in_maxgap = a gene desert; the scan only samples the ends\n")

## ---- per-lineage runs test (now correctly attributed) ----
runs_z <- function(s) {
  s <- s[!is.na(s)]; n1 <- sum(s == "A"); n2 <- sum(s == "B"); n <- n1 + n2
  if (n1 < 5 || n2 < 5) return(NA_real_)
  r <- 1 + sum(s[-1] != s[-n]); mu <- 2*n1*n2/n + 1
  vv <- 2*n1*n2*(2*n1*n2 - n)/(n^2*(n-1)); (r - mu)/sqrt(vv) }

lin <- v %>% arrange(pair, posA) %>% group_by(pair, lineage) %>%
  filter(n() >= 20) %>% ungroup()
per_lin <- lin %>% group_by(pair, lineage) %>%
  summarise(n = n(), frac_A = mean(vote == "A"), z = runs_z(vote), .groups = "drop") %>%
  mutate(p = 2*pnorm(-abs(z))) %>% arrange(z)
cat(sprintf("\n=== per-lineage runs test (%d series with >=20 votes) ===\n", nrow(per_lin)))
print(as.data.frame(head(per_lin, 12)), digits = 3)
cat(sprintf("... z summary: min %.2f  median %.2f  max %.2f | negative = clustered = tracts\n",
            min(per_lin$z, na.rm=TRUE), median(per_lin$z, na.rm=TRUE), max(per_lin$z, na.rm=TRUE)))

## ---- cross-lineage concordance vs an EXPLICIT null ----
gene <- v %>% group_by(pair, anchor, posA) %>%
  summarise(nlin = n_distinct(lineage), fA = mean(vote == "A"),
            agree = max(fA, 1-fA), call = ifelse(fA > .5, "A", ifelse(fA < .5, "B", NA)),
            .groups = "drop")
cat("\n=== lineages per gene ===\n"); print(as.data.frame(count(gene, nlin)))

null_agree <- function(nl, pA, B = 999) {
  mean(replicate(B, mean(sapply(nl, function(k) {
    x <- rbinom(1, k, pA); max(x, k-x)/k }))))
}
conc <- gene %>% filter(nlin >= 2) %>% group_by(pair) %>%
  summarise(n = n(), mean_agree = mean(agree), frac_unanimous = mean(agree == 1),
            pA = mean(fA), nl = list(nlin), .groups = "drop") %>%
  rowwise() %>%
  mutate(null_mean_agree = null_agree(nl[[1]], pA, B),
         excess = mean_agree - null_mean_agree) %>% ungroup() %>% select(-nl)
cat("\n=== cross-lineage concordance vs simulated null ===\n")
print(as.data.frame(conc), digits = 3)
cat("excess >> 0 => lineages move together => flips are DIONAEA breakpoints\n")
cat("excess ~ 0  => votes independent => no shared tract structure detectable\n")

## ---- plot, masking windows that span a gene desert ----
g <- gene %>% filter(nlin >= 1, !is.na(call)) %>% arrange(pair, posA)
sw <- g %>% group_by(pair) %>% group_modify(function(d, k) {
  n <- nrow(d); if (n < WIN) return(tibble())
  s <- seq_len(n - WIN + 1)
  tibble(mid_Mb = sapply(s, function(i) median(d$posA[i:(i+WIN-1)]))/1e6,
         span   = sapply(s, function(i) d$posA[i+WIN-1] - d$posA[i]),
         fA     = sapply(s, function(i) mean(d$call[i:(i+WIN-1)] == "A"))) }) %>%
  ungroup() %>% filter(span <= MAXGAP)

p <- ggplot(sw, aes(mid_Mb, fA)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
  geom_line(colour = "#7570b3", linewidth = .7) +
  geom_point(data = g, aes(posA/1e6, ifelse(call == "A", 1.06, -0.06)),
             shape = 124, size = 1.6, alpha = .5, colour = "#333333") +
  facet_wrap(~ pair, ncol = 2, scales = "free_x") +
  coord_cartesian(ylim = c(-0.1, 1.1)) +
  labs(title = "Which Dionaea homeolog do Drosera copies pair with?",
       subtitle = sprintf("%d-gene sliding mean, windows spanning >%.0f Mb suppressed; ticks = genes (top sgA, bottom sgB)",
                          WIN, MAXGAP/1e6),
       x = "position on sgA chromosome (Mb)", y = "fraction voting sgA") +
  theme_bw(base_size = 10)
ggsave("FIG13d_tract_scan7iq.png", p, width = 12, height = 11, dpi = 170, device = agg_png)
write_csv(gene, "tract_gene_calls7iq.csv"); write_csv(per_lin, "tract_per_lineage7iq.csv")
cat("\nWROTE: FIG13d_tract_scan7iq.png tract_votes7iq.csv tract_gene_calls7iq.csv tract_per_lineage7iq.csv\n")
