#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# One figure per Dionaea chromosome pair.
#   columns = position-ordered windows of consecutive anchor genes
#   rows    = Drosera lineage blocks, GROUPED and BORDER-COLOURED by species
#   fill    = fraction of that lineage's votes in that window going to sgA
#   vertical lines = GENESPACE Dionaea block boundaries (does a flip fall inside a block?)
#   bottom rows = per-species consensus, then overall consensus
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
W <- 8; MINROW <- 6

SPC <- c(regia = "#5B4EA8", binata = "#12795E", paradoxa = "#C0392B",
         scorpioides = "#B7791F", capensis = "#2C7FB8")
ABB <- c(regia = "reg", binata = "bin", paradoxa = "par",
         scorpioides = "sco", capensis = "cap")

blk <- read_csv(file.path(GSD, "syntenicBlock_coordinates.csv"), show_col_types = FALSE)
v <- read_csv("tract_votes_blocks6.csv", show_col_types = FALSE) %>%
  mutate(species = sub("Drosera_", "", sp),
         unit = sprintf("%s %s b%s", ABB[species],
                        sub("_hap1$|_collapsed$", "", sub("^chr", "c", lin_chr)), blk))

dnb <- blk %>% filter(genome1 == "Nepenthes_gracilis" | genome2 == "Nepenthes_gracilis") %>%
  mutate(n1 = genome1 == "Nepenthes_gracilis",
         spg = ifelse(n1, genome2, genome1), nep_chr = ifelse(n1, chr1, chr2),
         sp_chr = ifelse(n1, chr2, chr1),
         s = ifelse(n1, pmin(startBp2, endBp2), pmin(startBp1, endBp1)),
         e = ifelse(n1, pmax(startBp2, endBp2), pmax(startBp1, endBp1)),
         nhit = ifelse(n1, nHits2, nHits1)) %>%
  filter(spg == "Dionaea_muscipula", grepl("_dom$", nep_chr), !is.na(s)) %>%
  select(blkID, nep_chr, sp_chr, s, e, nhit)

dblk <- function(chr_, pos_, nchr_) {
  c <- dnb[dnb$sp_chr == chr_ & dnb$nep_chr == nchr_ & dnb$s <= pos_ & dnb$e >= pos_, ]
  if (!nrow(c)) return(NA_character_)
  as.character(c$blkID[which.max(c$nhit)])
}
anch <- v %>% distinct(pair, chrA, chrB, nep_chr, anchor, posA)
anch$dioblk <- mapply(dblk, anch$chrA, anch$posA, anch$nep_chr)
anch <- anch %>% arrange(pair, posA) %>% group_by(pair) %>%
  mutate(win = (row_number() - 1) %/% W + 1) %>% ungroup()

wsum <- anch %>% group_by(pair, chrA, chrB, nep_chr, win) %>%
  summarise(n_anch = n(), mid_Mb = median(posA)/1e6,
            dioblk = names(sort(table(dioblk), decreasing = TRUE))[1], .groups = "drop") %>%
  arrange(pair, win)

V <- v %>% left_join(select(anch, pair, anchor, win), by = c("pair", "anchor"))

M <- V %>% group_by(pair, win, species, unit) %>%
  summarise(n = n(), att = mean(vote == "A"), .groups = "drop") %>%
  group_by(pair, unit) %>% filter(sum(n) >= MINROW, n_distinct(win) >= 2) %>% ungroup()

## per-species and overall consensus rows
SC <- V %>% group_by(pair, win, species) %>%
  summarise(n = n(), att = mean(vote == "A"), .groups = "drop") %>%
  mutate(unit = paste0("== ", ABB[species], " consensus =="))
OC <- V %>% group_by(pair, win) %>%
  summarise(n = n(), att = mean(vote == "A"), .groups = "drop") %>%
  mutate(species = "ALL", unit = "=== ALL consensus ===")
ALL <- bind_rows(M, SC, OC)

cat(sprintf("rows: %d lineage blocks + %d consensus rows across %d pairs\n",
            n_distinct(paste(M$pair, M$unit)), n_distinct(paste(SC$pair, SC$unit)) + n_distinct(OC$pair),
            n_distinct(ALL$pair)))
cat("\nlineage blocks per species per Dionaea pair (copy structure):\n")
print(as.data.frame(M %>% distinct(pair, species, unit) %>% count(pair, species) %>%
                      pivot_wider(names_from = species, values_from = n, values_fill = 0)))

dir.create("figs_window", showWarnings = FALSE)
for (pp in sort(unique(ALL$pair))) {
  d  <- filter(ALL, pair == pp)
  ws <- filter(wsum, pair == pp)
  ord <- d %>% filter(species != "ALL", !grepl("consensus", unit)) %>%
    group_by(species, unit) %>% summarise(m = mean(att), .groups = "drop") %>%
    mutate(species = factor(species, names(SPC))) %>%
    arrange(species, desc(m)) %>% pull(unit)
  cons <- c(paste0("== ", ABB[names(SPC)], " consensus =="), "=== ALL consensus ===")
  cons <- cons[cons %in% d$unit]
  d$unit <- factor(d$unit, levels = rev(c(ord, cons)))
  d$species <- factor(d$species, c(names(SPC), "ALL"))

  bb <- ws %>% mutate(chg = dioblk != lag(dioblk)) %>% filter(chg) %>% pull(win)
  lab <- setNames(sprintf("%d\n%.0f", ws$win, ws$mid_Mb), ws$win)

  ttl <- sprintf("Nepenthes %s   |   Dionaea %s (sgA)  vs  %s (sgB)",
                 ws$nep_chr[1], ws$chrA[1], ws$chrB[1])
  nr <- n_distinct(d$unit)
  p <- ggplot(d, aes(factor(win), unit, fill = att, colour = species)) +
    geom_tile(linewidth = .7, width = .95, height = .88) +
    geom_text(aes(label = n), size = 2.1, colour = "grey12") +
    { if (length(bb)) geom_vline(xintercept = bb - 0.5, colour = "black",
                                 linetype = 2, linewidth = .5) } +
    geom_hline(yintercept = length(cons) + 0.5, colour = "black", linewidth = .6) +
    scale_fill_gradient2(low = "#d95f02", mid = "grey94", high = "#1b9e77",
                         midpoint = .5, limits = c(0, 1), na.value = "white") +
    scale_colour_manual(values = c(SPC, ALL = "black"), drop = FALSE) +
    scale_x_discrete(labels = lab) +
    labs(title = ttl,
         subtitle = sprintf("columns = %d-anchor windows left to right (window no. / Mb); dashed lines = GENESPACE Dionaea block boundaries\nrow border colour = species; number in cell = votes; bottom rows = consensus", W),
         x = "window along Dionaea sgA", y = NULL,
         fill = "frac sgA", colour = "species") +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(size = 6.5),
          axis.text.x = element_text(size = 6),
          panel.grid = element_blank())
  ggsave(sprintf("figs_window/FIG25_%s.png", pp), p,
         width = max(8, 0.42 * n_distinct(d$win) + 5),
         height = max(5, 0.19 * nr + 2.4), dpi = 175, limitsize = FALSE, device = agg_png)
}
write_csv(ALL, "window_matrix_rows.csv"); write_csv(wsum, "window_index2.csv")
cat(sprintf("\nWROTE: figs_window/FIG25_<pair>.png (%d files)\n", n_distinct(ALL$pair)))
