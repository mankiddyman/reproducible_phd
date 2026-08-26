#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# FIG17 but with Dionaea tiled into POSITION-ORDERED windows of consecutive anchors
# instead of GENESPACE blocks. Rows = Drosera lineage-blocks. Cell = fraction of that
# lineage's votes in that window going to sgA.
# A colour flip between adjacent columns = candidate breakpoint.
# Column strip shows which GENESPACE Dionaea block each window sits in, so you can see
# whether a break falls INSIDE a block (block must be split) or at its boundary.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ragg); library(patchwork)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
W <- 8; MINROW <- 6

blk <- read_csv(file.path(GSD, "syntenicBlock_coordinates.csv"), show_col_types = FALSE)
v <- read_csv("tract_votes_blocks6.csv", show_col_types = FALSE) %>%
  mutate(unit = paste0(sub("Drosera_", "", sp) %>% substr(1,3), " ",
                       sub("_hap1$|_collapsed$", "", sub("^chr","c",lin_chr)), " b", blk))

## Dionaea blocks (sgA side) vs Nepenthes
dnb <- blk %>% filter(genome1 == "Nepenthes_gracilis" | genome2 == "Nepenthes_gracilis") %>%
  mutate(n1 = genome1 == "Nepenthes_gracilis",
         sp = ifelse(n1, genome2, genome1), nep_chr = ifelse(n1, chr1, chr2),
         sp_chr = ifelse(n1, chr2, chr1),
         s = ifelse(n1, pmin(startBp2, endBp2), pmin(startBp1, endBp1)),
         e = ifelse(n1, pmax(startBp2, endBp2), pmax(startBp1, endBp1)),
         nhit = ifelse(n1, nHits2, nHits1)) %>%
  filter(sp == "Dionaea_muscipula", grepl("_dom$", nep_chr), !is.na(s)) %>%
  select(blkID, nep_chr, sp_chr, s, e, nhit)

dblk <- function(chr_, pos_, nchr_) {
  c <- dnb[dnb$sp_chr == chr_ & dnb$nep_chr == nchr_ & dnb$s <= pos_ & dnb$e >= pos_, ]
  if (!nrow(c)) return(NA_character_)
  as.character(c$blkID[which.max(c$nhit)])
}
anch <- v %>% distinct(pair, chrA, nep_chr, anchor, posA)
anch$dioblk <- mapply(dblk, anch$chrA, anch$posA, anch$nep_chr)
cat(sprintf("anchors: %d | with a Dionaea block: %.1f%%\n",
            nrow(anch), 100*mean(!is.na(anch$dioblk))))

## tile into windows of W consecutive anchors, ordered by position
anch <- anch %>% arrange(pair, posA) %>% group_by(pair) %>%
  mutate(idx = row_number(), win = (idx - 1) %/% W + 1) %>% ungroup()
wsum <- anch %>% group_by(pair, win) %>%
  summarise(n_anch = n(), start_Mb = min(posA)/1e6, end_Mb = max(posA)/1e6,
            mid_Mb = median(posA)/1e6,
            dioblks = paste(sort(unique(na.omit(dioblk))), collapse = ","),
            .groups = "drop")
cat(sprintf("windows (W=%d anchors): %d across %d pairs\n", W, nrow(wsum), n_distinct(wsum$pair)))

V <- v %>% left_join(select(anch, pair, anchor, win), by = c("pair","anchor"))

M <- V %>% group_by(pair, win, unit) %>%
  summarise(n = n(), att = mean(vote == "A"), .groups = "drop") %>%
  group_by(pair, unit) %>% filter(sum(n) >= MINROW, n_distinct(win) >= 2) %>% ungroup() %>%
  left_join(select(wsum, pair, win, mid_Mb, dioblks), by = c("pair","win"))
cat(sprintf("matrix cells: %d | rows kept: %d\n", nrow(M), n_distinct(paste(M$pair, M$unit))))

## ---- is each Dionaea GENESPACE block internally consistent? ----------------
cons <- V %>% left_join(select(anch, pair, anchor, dioblk), by = c("pair","anchor")) %>%
  filter(!is.na(dioblk)) %>%
  group_by(pair, dioblk, win) %>%
  summarise(n = n(), att = mean(vote == "A"), .groups = "drop") %>%
  group_by(pair, dioblk) %>% filter(n_distinct(win) >= 2, sum(n) >= 15) %>%
  summarise(n_win = n_distinct(win), tot = sum(n),
            min_att = min(att), max_att = max(att), spread = max_att - min_att,
            crosses = min_att < 0.4 & max_att > 0.6, .groups = "drop") %>%
  arrange(desc(spread))
cat("\n=== Dionaea GENESPACE blocks spanning >=2 windows: internally consistent? ===\n")
print(as.data.frame(cons), digits = 3)
if (nrow(cons)) cat(sprintf("\n%d/%d blocks contain windows on OPPOSITE sides => candidate split points\n",
                            sum(cons$crosses), nrow(cons)))
write_csv(cons, "dioblock_internal_consistency.csv")
write_csv(M, "window_matrix.csv"); write_csv(wsum, "window_index.csv")

## ---- plot ------------------------------------------------------------------
lab <- wsum %>% mutate(collab = sprintf("%02d\n%.0f", win, mid_Mb))
M <- M %>% left_join(select(lab, pair, win, collab), by = c("pair","win"))
M$collab <- factor(M$collab, levels = unique(lab$collab[order(lab$pair, lab$win)]))

p1 <- ggplot(M, aes(reorder(collab, win), reorder(unit, att), fill = att)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = n), size = 2, colour = "grey15") +
  facet_wrap(~ pair, scales = "free", ncol = 2) +
  scale_fill_gradient2(low = "#d95f02", mid = "grey93", high = "#1b9e77",
                       midpoint = .5, limits = c(0,1)) +
  labs(title = "Dionaea tiled into position-ordered windows of consecutive anchor genes",
       subtitle = sprintf("columns = windows of %d anchors, left to right along the chromosome (label: window no. / Mb)\nrows = Drosera lineage blocks; a colour flip across adjacent columns is a candidate breakpoint", W),
       x = "window along Dionaea sgA", y = "Drosera lineage block", fill = "frac sgA") +
  theme_bw(base_size = 7) +
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5))
ggsave("FIG23_window_matrix.png", p1, width = 15, height = 16, dpi = 165, device = agg_png)

## column strip: which GENESPACE block each window belongs to
p2 <- wsum %>% mutate(blkgrp = sub(",.*", "", dioblks)) %>%
  ggplot(aes(factor(win), 1, fill = blkgrp)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = blkgrp), size = 1.8, angle = 90) +
  facet_wrap(~ pair, scales = "free_x", ncol = 2) +
  labs(title = "Which GENESPACE Dionaea block each window falls in",
       subtitle = "a breakpoint inside one of these = that block must be split",
       x = "window", y = NULL) +
  theme_bw(base_size = 7) +
  theme(legend.position = "none", axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.text.x = element_text(size = 4))
ggsave("FIG24_window_blocks.png", p2, width = 13, height = 7, dpi = 165, device = agg_png)
cat("\nWROTE: FIG23_window_matrix.png FIG24_window_blocks.png window_matrix.csv dioblock_internal_consistency.csv\n")
