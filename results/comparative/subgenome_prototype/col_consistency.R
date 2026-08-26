#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Do different Dionaea block pairs along one chromosome pair give the SAME partition
# of Drosera blocks? Same => no HE, that stretch is one subgenome.
# Inverted => the two Dionaea chromosomes swapped ancestry between those positions.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
MINSHARE <- 2   # Drosera blocks in common between two columns

A <- read_csv("block_association.csv", show_col_types = FALSE) %>%
  mutate(DA = sub("Dionaea_muscipula_vs_Nepenthes_gracilis: ", "", DA),
         DB = sub("Dionaea_muscipula_vs_Nepenthes_gracilis: ", "", DB),
         dros = sub("Drosera_[a-z]+_vs_Nepenthes_gracilis: ", "", dros),
         col = paste(DA, DB, sep = "|"))

## column position = midpoint of the DA block on its chromosome
blk <- read_csv(file.path(GSD, "syntenicBlock_coordinates.csv"), show_col_types = FALSE)
nbd <- blk %>% filter(genome1 == "Nepenthes_gracilis" | genome2 == "Nepenthes_gracilis") %>%
  mutate(n1 = genome1 == "Nepenthes_gracilis",
         sp = ifelse(n1, genome2, genome1),
         sp_chr = ifelse(n1, chr2, chr1),
         s = ifelse(n1, pmin(startBp2, endBp2), pmin(startBp1, endBp1)),
         e = ifelse(n1, pmax(startBp2, endBp2), pmax(startBp1, endBp1))) %>%
  filter(sp == "Dionaea_muscipula") %>%
  transmute(bid = as.character(blkID),
            key = paste0(sub("_sg[0-9]+_s", "-s", sub("^chr","c",sp_chr)), "#", blkID),
            mid = (s + e)/2, s, e)
A <- A %>% left_join(nbd %>% transmute(DA = key, posA = mid, sA = s, eA = e), by = "DA")

# use ALL cells with >=3 votes; keep the signed deviation so weak cells count less
conf <- A %>% filter(n >= 3) %>% mutate(dev = att - 0.5, side = sign(dev)) %>%
  filter(side != 0)
cat(sprintf("cells used: %d over %d columns, %d Drosera blocks (mean %.1f cells/column)\n",
            nrow(conf), n_distinct(conf$col), n_distinct(conf$dros),
            nrow(conf)/n_distinct(conf$col)))
sh <- conf %>% group_by(nep_chr) %>%
  summarise(cols = n_distinct(col), cells = n(), .groups = "drop")
cat("\ncells per region:\n"); print(as.data.frame(sh))

## ---- pairwise column comparison within each Nepenthes region ---------------
res <- bind_rows(lapply(unique(conf$nep_chr), function(rg) {
  d <- filter(conf, nep_chr == rg)
  cols <- sort(unique(d$col)); if (length(cols) < 2) return(NULL)
  cb <- combn(cols, 2)
  bind_rows(lapply(seq_len(ncol(cb)), function(i) {
    x <- filter(d, col == cb[1,i]); y <- filter(d, col == cb[2,i])
    sh <- intersect(x$dros, y$dros); if (length(sh) < MINSHARE) return(NULL)
    sx <- x$side[match(sh, x$dros)]; sy <- y$side[match(sh, y$dros)]
    wx <- abs(x$dev[match(sh, x$dros)]); wy <- abs(y$dev[match(sh, y$dros)])
    w  <- wx * wy
    fa <- if (sum(w) > 0) sum(w * (sx == sy))/sum(w) else mean(sx == sy)
    tibble(nep_chr = rg, col1 = cb[1,i], col2 = cb[2,i],
           pos1 = x$posA[1]/1e6, pos2 = y$posA[1]/1e6,
           n_shared = length(sh), n_agree = sum(sx == sy),
           frac_agree = mean(sx == sy), w_agree = fa,
           verdict = ifelse(fa >= .8, "SAME", ifelse(fa <= .2, "INVERTED", "mixed")))
  }))
}))
if (!nrow(res)) {
  cat(sprintf("\nNO comparable column pairs at MINSHARE=%d.\n", MINSHARE))
  cat("Columns do not share enough Drosera blocks -- the block grid is too fine.\n")
  ov <- conf %>% group_by(nep_chr, dros) %>% summarise(ncol = n_distinct(col), .groups="drop")
  cat("\nDrosera blocks by how many Dionaea columns they appear in:\n")
  print(as.data.frame(count(ov, ncol)))
  quit(save = "no")
}
write_csv(res, "column_consistency.csv")

cat(sprintf("\ncomparable column pairs (>=%d shared Drosera blocks): %d\n", MINSHARE, nrow(res)))
cat("\n=== verdicts ===\n"); print(as.data.frame(count(res, verdict)))
cat("\nSAME = no HE between those two positions | INVERTED = HE breakpoint between them\n")

cat("\n=== all comparisons, by region, ordered by position ===\n")
for (rg in sort(unique(res$nep_chr))) {
  cat(sprintf("\n-- %s --\n", rg))
  print(as.data.frame(filter(res, nep_chr == rg) %>%
    arrange(pos1, pos2) %>%
    transmute(col1, pos1_Mb = round(pos1,1), col2, pos2_Mb = round(pos2,1),
              n_shared, frac_agree = round(frac_agree,2),
              w_agree = round(w_agree,2), verdict)), right = FALSE)
}

inv <- filter(res, verdict == "INVERTED")
if (nrow(inv)) {
  cat("\n=== INVERTED column pairs — candidate HE breakpoints ===\n")
  print(as.data.frame(inv %>% transmute(nep_chr, col1, pos1_Mb = round(pos1,1),
                                        col2, pos2_Mb = round(pos2,1),
                                        n_shared, frac_agree = round(frac_agree,2))),
        right = FALSE)
} else cat("\nno INVERTED column pairs: no HE detected between Dionaea homeologs at this resolution\n")

p <- ggplot(res, aes(w_agree, fill = verdict)) +
  geom_histogram(binwidth = 0.1, colour = "white") +
  geom_vline(xintercept = c(0.2, 0.8), linetype = 2, colour = "grey40") +
  scale_fill_manual(values = c(SAME = "#1b9e77", INVERTED = "#d95f02", mixed = "grey70")) +
  labs(title = "Do two Dionaea block pairs give the same partition of Drosera blocks?",
       subtitle = "1 = identical partition (no HE between them); 0 = inverted (HE breakpoint)",
       x = "fraction of shared Drosera blocks agreeing", y = "column pairs") +
  theme_bw(base_size = 11)
ggsave("FIG18_column_consistency.png", p, width = 8, height = 5, dpi = 180, device = agg_png)
cat("\nWROTE: FIG18_column_consistency.png column_consistency.csv\n")
