#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# One multi-page PDF, a page per Dionaea chromosome pair.
#   columns = position-ordered windows of consecutive anchor genes
#   rows    = Drosera lineage blocks, grouped by species with horizontal dividers
#   cell    = STACKED PROPORTION: green = share of votes for sgA, orange = sgB
#   dashed vertical lines = GENESPACE Dionaea block boundaries
#   bottom band = per-species consensus, then ALL consensus
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
W <- 8; MINROW <- 6
SPORD <- c("regia","binata","paradoxa","scorpioides","capensis")
ABB <- c(regia="reg", binata="bin", paradoxa="par", scorpioides="sco", capensis="cap")

blk <- read_csv(file.path(GSD,"syntenicBlock_coordinates.csv"), show_col_types = FALSE)
v <- read_csv("tract_votes_blocks6.csv", show_col_types = FALSE) %>%
  mutate(species = sub("Drosera_","",sp),
         unit = sprintf("%s %s b%s", ABB[species],
                        sub("_hap1$|_collapsed$","", sub("^chr","c",lin_chr)), blk))

dnb <- blk %>% filter(genome1=="Nepenthes_gracilis" | genome2=="Nepenthes_gracilis") %>%
  mutate(n1 = genome1=="Nepenthes_gracilis",
         spg = ifelse(n1,genome2,genome1), nep_chr = ifelse(n1,chr1,chr2),
         sp_chr = ifelse(n1,chr2,chr1),
         s = ifelse(n1,pmin(startBp2,endBp2),pmin(startBp1,endBp1)),
         e = ifelse(n1,pmax(startBp2,endBp2),pmax(startBp1,endBp1)),
         nhit = ifelse(n1,nHits2,nHits1)) %>%
  filter(spg=="Dionaea_muscipula", grepl("_dom$",nep_chr), !is.na(s)) %>%
  select(blkID, nep_chr, sp_chr, s, e, nhit)
dbf <- function(c_,p_,n_) { r <- dnb[dnb$sp_chr==c_ & dnb$nep_chr==n_ & dnb$s<=p_ & dnb$e>=p_,]
  if(!nrow(r)) NA_character_ else as.character(r$blkID[which.max(r$nhit)]) }

anch <- v %>% distinct(pair, chrA, chrB, nep_chr, anchor, posA)
anch$dioblk <- mapply(dbf, anch$chrA, anch$posA, anch$nep_chr)
anch <- anch %>% arrange(pair, posA) %>% group_by(pair) %>%
  mutate(win = (row_number()-1) %/% W + 1) %>% ungroup()
wsum <- anch %>% group_by(pair, chrA, chrB, nep_chr, win) %>%
  summarise(mid_Mb = median(posA)/1e6,
            dioblk = names(sort(table(dioblk), decreasing=TRUE))[1], .groups="drop") %>%
  arrange(pair, win)

V <- v %>% left_join(select(anch, pair, anchor, win), by = c("pair","anchor"))
M  <- V %>% count(pair, win, species, unit, vote, name="k") %>%
  group_by(pair, unit) %>% filter(sum(k) >= MINROW, n_distinct(win) >= 2) %>% ungroup()
SC <- V %>% count(pair, win, species, vote, name="k") %>%
  mutate(unit = paste0("~ ", ABB[species], " consensus"))
OC <- V %>% count(pair, win, vote, name="k") %>%
  mutate(species = "consensus", unit = "~~ ALL consensus")
D <- bind_rows(M, SC, OC) %>%
  group_by(pair, win, unit, species) %>%
  mutate(frac = k/sum(k), ntot = sum(k)) %>% ungroup() %>%
  mutate(vote = factor(vote, c("B","A")))

pdf("FIG26_window_matrix_all.pdf", width = 16, height = 11, onefile = TRUE)
for (pp in sort(unique(D$pair))) {
  d <- filter(D, pair == pp); ws <- filter(wsum, pair == pp)

  rows <- d %>% distinct(species, unit) %>%
    mutate(grp = ifelse(species == "consensus", "consensus", species)) %>%
    left_join(d %>% filter(vote=="A") %>% group_by(unit) %>%
                summarise(m = weighted.mean(frac, ntot), .groups="drop"), by="unit") %>%
    mutate(m = replace_na(m, .5),
           iscons = grepl("consensus", unit),
           spf = factor(ifelse(iscons, "zz_consensus", species),
                        c(SPORD, "zz_consensus"))) %>%
    arrange(spf, desc(m))
  d$unit <- factor(d$unit, levels = rev(rows$unit))

  # y positions of species-group dividers
  rr <- rows %>% mutate(y = rev(seq_len(n())))
  div <- rr %>% group_by(spf) %>% summarise(top = max(y), .groups="drop") %>%
    arrange(desc(top)) %>% mutate(cut = top + 0.5) %>% pull(cut)
  div <- div[div < nrow(rr)]
  labs_sp <- rr %>% group_by(spf) %>%
    summarise(y = mean(range(y)), lab = ifelse(spf[1]=="zz_consensus","CONSENSUS",
                                               as.character(spf[1])), .groups="drop")

  bb <- ws %>% mutate(chg = dioblk != lag(dioblk)) %>% filter(chg) %>% pull(win)
  xlab <- setNames(sprintf("%d\n%.0f", ws$win, ws$mid_Mb), ws$win)

  # ggplot cannot stack inside a heatmap cell; draw with geom_rect instead
  dd <- d %>% arrange(pair, unit, win, vote) %>%
    group_by(pair, unit, win) %>%
    mutate(ymaxf = cumsum(frac), yminf = ymaxf - frac,
           yi = as.integer(unit),
           ymin = yi - 0.45 + 0.9*yminf, ymax = yi - 0.45 + 0.9*ymaxf,
           xi = match(win, ws$win), xmin = xi - 0.46, xmax = xi + 0.46) %>% ungroup()

  p <- ggplot(dd) +
    geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=vote), colour=NA) +
    geom_rect(data = distinct(dd, unit, win, xmin, xmax, yi, ntot),
              aes(xmin=xmin, xmax=xmax, ymin=yi-0.45, ymax=yi+0.45),
              fill=NA, colour="grey55", linewidth=.25) +
    geom_text(data = distinct(dd, unit, win, xi, yi, ntot),
              aes(xi, yi, label = ntot), size = 2, colour = "grey15") +
    { if (length(bb)) geom_vline(xintercept = match(bb, ws$win) - .5,
                                 linetype = 2, colour = "black", linewidth = .5) } +
    geom_hline(yintercept = div, colour = "black", linewidth = .7) +
    geom_text(data = labs_sp, aes(x = -1.6, y = y, label = lab),
              angle = 90, size = 3.2, fontface = 2, inherit.aes = FALSE) +
    scale_fill_manual(values = c(A = "#1b9e77", B = "#d95f02"),
                      labels = c(A = "sgA", B = "sgB"), breaks = c("A","B")) +
    scale_y_continuous(breaks = seq_len(nrow(rr)),
                       labels = levels(d$unit), expand = expansion(add = .6)) +
    scale_x_continuous(breaks = seq_along(ws$win), labels = unname(xlab[as.character(ws$win)]),
                       limits = c(-2.6, length(ws$win) + .6), expand = c(0,0)) +
    labs(title = sprintf("Nepenthes %s   |   Dionaea %s (sgA) vs %s (sgB)",
                         ws$nep_chr[1], ws$chrA[1], ws$chrB[1]),
         subtitle = sprintf("each cell is a stacked proportion of that block's votes in that window: green = sgA, orange = sgB (number = total votes)\ncolumns = %d-anchor windows left to right (window no. / Mb); dashed lines = GENESPACE Dionaea block boundaries", W),
         x = "window along Dionaea sgA", y = NULL, fill = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 6.2),
          axis.text.x = element_text(size = 6),
          legend.position = "top")
  print(p)
}
dev.off()
cat("WROTE: FIG26_window_matrix_all.pdf (one page per Dionaea chromosome pair)\n")
