#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# BLOCK-level allyship. A Drosera chromosome can carry several syntenic blocks against
# one Nepenthes chromosome, and those blocks may descend from different Dionaea
# subgenomes. So the unit is (species, chromosome, syntenic block), not chromosome.
# Key test: do blocks on the SAME chromosome disagree? If yes, chromosome-level calls
# are invalid. If no, they were fine.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
MINN <- 8

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
blk <- read_csv(file.path(GSD, "syntenicBlock_coordinates.csv"), show_col_types = FALSE)
anch <- read_csv("synteny_ortho_hits.csv", show_col_types = FALSE) %>% distinct(nep_gene, nep_chr)
v <- read_csv("tract_votes7.csv", show_col_types = FALSE)

## Nepenthes-anchored blocks, species-side coordinates
nb <- blk %>% filter(genome1 == "Nepenthes_gracilis" | genome2 == "Nepenthes_gracilis") %>%
  mutate(n1 = genome1 == "Nepenthes_gracilis",
         sp = ifelse(n1, genome2, genome1), nep_chr = ifelse(n1, chr1, chr2),
         sp_chr = ifelse(n1, chr2, chr1),
         sp_s = ifelse(n1, pmin(startBp2, endBp2), pmin(startBp1, endBp1)),
         sp_e = ifelse(n1, pmax(startBp2, endBp2), pmax(startBp1, endBp1)),
         nhit = ifelse(n1, nHits2, nHits1)) %>%
  filter(sp != "Nepenthes_gracilis", grepl("_dom$", nep_chr), !is.na(sp_s)) %>%
  select(blkID, sp, nep_chr, sp_chr, sp_s, sp_e, nhit)
cat(sprintf("Nepenthes-anchored blocks: %d\n", nrow(nb)))

## voting gene -> block
gp <- bed %>% transmute(sp_gene = id, genome, gmid = (start + end)/2)
v <- v %>% left_join(gp, by = c("sp_gene", "sp" = "genome")) %>%
  left_join(anch, by = c("anchor" = "nep_gene"))

assign_blk <- function(sp_, chr_, mid_, nchr_) {
  if (is.na(mid_) || is.na(nchr_)) return(NA_character_)
  c <- nb[nb$sp == sp_ & nb$sp_chr == chr_ & nb$nep_chr == nchr_ &
          nb$sp_s <= mid_ & nb$sp_e >= mid_, ]
  if (!nrow(c)) return(NA_character_)
  as.character(c$blkID[which.max(c$nhit)])
}
v$blk <- mapply(assign_blk, v$sp, v$lin_chr, v$gmid, v$nep_chr)
cat(sprintf("votes: %d | assigned to a block: %.1f%%\n",
            nrow(v), 100*mean(!is.na(v$blk))))
v <- filter(v, !is.na(blk)) %>%
  mutate(unit = paste0(sub("Drosera_", "", sp), " ",
                       sub("_hap1$|_collapsed$", "", lin_chr), " b", blk))
write_csv(v, "tract_votes_blocks7.csv")

## ---- per-block skew ---------------------------------------------------------
B <- v %>% group_by(pair, chrA, chrB, sp, lin_chr, blk, unit) %>%
  summarise(n = n(), A = sum(vote == "A"), frac_A = A/n, .groups = "drop")
cat("\n=== votes per block ===\n")
print(summary(B$n))
cat(sprintf("blocks with >=%d votes: %d of %d\n", MINN, sum(B$n >= MINN), nrow(B)))

Bs <- B %>% filter(n >= MINN) %>% rowwise() %>%
  mutate(p = binom.test(A, n, 0.5)$p.value,
         lo = binom.test(A, n, 0.5)$conf.int[1],
         hi = binom.test(A, n, 0.5)$conf.int[2]) %>% ungroup() %>%
  mutate(p_adj = p.adjust(p, "BH"),
         call = case_when(p_adj < .05 & frac_A > .5 ~ "A",
                          p_adj < .05 & frac_A < .5 ~ "B", TRUE ~ "?")) %>%
  arrange(pair, desc(frac_A))
write_csv(Bs, "block_skew7.csv")
cat(sprintf("\ncalled: %d A, %d B, %d ambiguous (of %d blocks)\n",
            sum(Bs$call=="A"), sum(Bs$call=="B"), sum(Bs$call=="?"), nrow(Bs)))

for (pp in unique(Bs$pair)) {
  d <- Bs[Bs$pair == pp, ]
  cat(sprintf("\n-- Dionaea %s  (A = %s / B = %s) --\n", pp, d$chrA[1], d$chrB[1]))
  print(as.data.frame(d[, c("unit","n","frac_A","lo","hi","p_adj","call")]), digits = 3)
}

## ---- THE TEST: do blocks on the same chromosome disagree? -------------------
multi <- B %>% filter(n >= 5) %>%
  group_by(pair, sp, lin_chr) %>% filter(n_distinct(blk) >= 2) %>%
  summarise(n_blk = n_distinct(blk), tot = sum(n),
            min_f = min(frac_A), max_f = max(frac_A), spread = max_f - min_f,
            crosses_half = min_f < .5 & max_f > .5, .groups = "drop") %>%
  arrange(desc(spread))
cat("\n=== same chromosome, >=2 blocks with >=5 votes each ===\n")
print(as.data.frame(multi), digits = 3)
if (nrow(multi)) cat(sprintf(
  "\n%d/%d chromosomes have blocks on OPPOSITE sides of 0.5; median spread %.3f\n",
  sum(multi$crosses_half), nrow(multi), median(multi$spread)))
cat("many crossing => chromosome-level calls are invalid, blocks are mandatory\n")
cat("few crossing  => blocks within a chromosome agree; chromosome level was adequate\n")

p1 <- ggplot(Bs, aes(frac_A, reorder(unit, frac_A), colour = call)) +
  geom_vline(xintercept = .5, linetype = 2, colour = "grey40") +
  geom_linerange(aes(xmin = lo, xmax = hi)) + geom_point(aes(size = n)) +
  facet_wrap(~ pair, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c(A = "#1b9e77", B = "#d95f02", `?` = "grey65")) +
  scale_size_continuous(range = c(1.2, 3)) +
  labs(title = "Which Dionaea homeolog does each Drosera SYNTENIC BLOCK ally with?",
       subtitle = "one row per (species, chromosome, block); 95% CI",
       x = "fraction of genes voting sgA", y = NULL, colour = "call", size = "genes") +
  theme_bw(base_size = 8)
ggsave("FIG15c_block_skew6.png", p1, width = 12, height = 13, dpi = 170, device = agg_png)
cat("\nWROTE: FIG15c_block_skew6.png block_skew7.csv tract_votes_blocks7.csv\n")
