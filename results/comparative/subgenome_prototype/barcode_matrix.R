#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# BARCODE version. Each cell = one window x one (species, Drosera chromosome).
# Inside the cell, every individual vote is drawn as its own vertical stripe,
# ordered LEFT-TO-RIGHT by position along the Dionaea sgA chromosome.
#   stripe colour  = which Dionaea homeolog that gene's tree put it with
#   number of stripes = number of votes (no separate count needed)
#   thin strip above = Drosera GENESPACE block identity, alternating grey shades,
#                      so a colour flip inside one grey band = that block is not a
#                      unit of single ancestry and needs splitting.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
W <- 8
SPORD <- c("regia","binata","paradoxa","scorpioides","capensis")
ABB <- c(regia="reg", binata="bin", paradoxa="par", scorpioides="sco", capensis="cap")
GREY <- c("#c8c8c8", "#8a8a8a")     # alternating Drosera block shades

blk <- read_csv(file.path(GSD,"syntenicBlock_coordinates.csv"), show_col_types = FALSE)
v <- read_csv("tract_votes_blocks7.csv", show_col_types = FALSE) %>%
  mutate(species = sub("Drosera_","",sp),
         bid  = sub(".*:\\s*","", blk),
         dchr = sub("_hap1$|_collapsed$","", sub("^chr","c", lin_chr)),
         unit = paste(ABB[species], dchr))
cat(sprintf("votes: %d | rows (species x Drosera chromosome): %d\n",
            nrow(v), n_distinct(paste(v$pair, v$unit))))

## Dionaea block per anchor, for the dashed boundaries
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
  if(!nrow(r)) NA_character_ else sub(".*:\\s*","",as.character(r$blkID[which.max(r$nhit)])) }

anch <- v %>% distinct(pair, chrA, chrB, nep_chr, anchor, posA)
anch$dioblk <- mapply(dbf, anch$chrA, anch$posA, anch$nep_chr)
anch <- anch %>% arrange(pair, posA) %>% group_by(pair) %>%
  mutate(win = (row_number()-1) %/% W + 1) %>% ungroup()
wsum <- anch %>% group_by(pair, chrA, chrB, nep_chr, win) %>%
  summarise(mid_Mb = median(posA)/1e6,
            dioblk = names(sort(table(dioblk), decreasing=TRUE))[1], .groups="drop") %>%
  arrange(pair, win)

V <- v %>% left_join(select(anch, pair, anchor, win), by = c("pair","anchor"))

## rows: real chromosomes, then per-species consensus, then ALL
ROW <- bind_rows(
  V %>% transmute(pair, win, species, unit, posA, vote, bid, kind = "chr"),
  V %>% transmute(pair, win, species, unit = paste0("~ ", ABB[species], " consensus"),
                  posA, vote, bid = NA_character_, kind = "spcons"),
  V %>% transmute(pair, win, species = "consensus", unit = "~~ ALL votes",
                  posA, vote, bid = NA_character_, kind = "allcons"))

dir.create("figs_barcode", showWarnings = FALSE)
pdf("FIG30_barcode_matrix.pdf", width = 20, height = 12, onefile = TRUE)
np <- 0
for (pp in sort(unique(ROW$pair))) {
  d  <- filter(ROW, pair == pp)
  ws <- filter(wsum, pair == pp)
  if (!nrow(d)) next

  ## row order: species groups, then consensus block at the bottom of the panel
  meta <- d %>% distinct(species, unit, kind) %>%
    left_join(d %>% group_by(unit) %>%
                summarise(fA = mean(vote == "A"), .groups="drop"), by = "unit") %>%
    mutate(spf = factor(ifelse(kind == "chr", species, "zz_consensus"),
                        c(SPORD, "zz_consensus"))) %>%
    filter(!is.na(spf)) %>% arrange(spf, kind, desc(fA))
  d <- filter(d, unit %in% meta$unit)
  d$unit <- factor(d$unit, levels = rev(meta$unit))

  ## stripe geometry: order votes by physical position inside each cell
  dd <- d %>% arrange(unit, win, posA, vote) %>%
    group_by(unit, win) %>%
    mutate(j = row_number(), n = n()) %>% ungroup() %>%
    mutate(yi = as.integer(unit), xi = match(win, ws$win)) %>%
    filter(!is.na(xi)) %>%
    mutate(x0 = xi - 0.46 + (j - 1) * 0.92 / n,
           x1 = xi - 0.46 + j * 0.92 / n,
           yb0 = yi - 0.45, yb1 = yi + 0.14,       # barcode band
           ys0 = yi + 0.17, ys1 = yi + 0.45)       # Drosera-block strip

  ## alternating grey per Drosera block, running along the whole row
  strip <- dd %>% filter(kind == "chr", !is.na(bid)) %>%
    arrange(unit, posA) %>% group_by(unit) %>%
    mutate(chg = bid != lag(bid), chg = ifelse(is.na(chg), TRUE, chg),
           gi = cumsum(chg), shade = GREY[(gi - 1) %% 2 + 1]) %>% ungroup()

  bb <- ws %>% mutate(chg = dioblk != lag(dioblk)) %>% filter(chg) %>% pull(win)
  rr <- meta %>% mutate(y = rev(seq_len(n())))
  div <- rr %>% group_by(spf) %>% summarise(top = max(y), .groups="drop") %>%
    mutate(cut = top + 0.5) %>% pull(cut)
  div <- div[div > 0.5 & div < nrow(rr)]
  lab_sp <- rr %>% group_by(spf) %>%
    summarise(y = mean(range(y)),
              lab = ifelse(spf[1]=="zz_consensus","CONSENSUS", as.character(spf[1])),
              .groups="drop")

  p <- ggplot() +
    geom_rect(data = strip, aes(xmin=x0, xmax=x1, ymin=ys0, ymax=ys1, fill=shade)) +
    scale_fill_identity() +
    ggnewscale::new_scale_fill() +
    geom_rect(data = dd, aes(xmin=x0, xmax=x1, ymin=yb0, ymax=yb1, fill=vote)) +
    geom_rect(data = distinct(dd, unit, win, xi, yi),
              aes(xmin=xi-0.46, xmax=xi+0.46, ymin=yi-0.45, ymax=yi+0.45),
              fill=NA, colour="grey35", linewidth=.22) +
    { if (length(bb)) geom_vline(xintercept = match(bb, ws$win) - .5,
                                 linetype = 2, colour = "black", linewidth = .55) } +
    geom_hline(yintercept = div, colour = "black", linewidth = .8) +
    geom_text(data = lab_sp, aes(x = -1.3, y = y, label = lab),
              angle = 90, size = 3.6, fontface = 2, inherit.aes = FALSE) +
    scale_fill_manual(values = c(A = "#1b9e77", B = "#d95f02"),
                      breaks = c("A","B"),
                      labels = c(A = sprintf("sgA  (%s)", ws$chrA[1]),
                                 B = sprintf("sgB  (%s)", ws$chrB[1])),
                      name = "gene tree places the Drosera copy with") +
    scale_y_continuous(breaks = seq_len(nrow(rr)), labels = levels(d$unit),
                       expand = expansion(add = .7)) +
    scale_x_continuous(breaks = seq_along(ws$win),
                       labels = sprintf("%d\n%.0f", ws$win, ws$mid_Mb),
                       limits = c(-2.2, length(ws$win) + .6), expand = c(0,0)) +
    labs(title = sprintf("Nepenthes %s   |   Dionaea %s (sgA)  vs  %s (sgB)",
                         ws$nep_chr[1], ws$chrA[1], ws$chrB[1]),
         subtitle = paste0(
           "Each stripe is ONE gene, ordered left to right by position on Dionaea sgA. Number of stripes = number of votes.\n",
           "Grey strip above each barcode = Drosera GENESPACE block; the shade alternates at every block boundary.\n",
           "A colour flip that lands INSIDE one grey band means that Drosera block is not a unit of single ancestry.\n",
           sprintf("Columns = windows of %d consecutive anchor genes (window no. / Mb). Dashed lines = Dionaea block boundaries.", W)),
         x = "window along Dionaea sgA", y = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 6.5),
          axis.text.x = element_text(size = 6),
          legend.position = "top")
  print(p); np <- np + 1
}
dev.off()
cat(sprintf("WROTE: FIG30_barcode_matrix.pdf (%d pages)\n", np))
