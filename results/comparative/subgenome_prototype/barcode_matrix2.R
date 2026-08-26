#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# GENE-ALIGNED barcode. Every row uses the same horizontal slot for the same anchor gene,
# so a slot can be read vertically across all species.
#   one SLOT   = one anchor gene (windows of W anchors, left to right along Dionaea sgA)
#   one BAND   = one gene copy of that species; green = placed with sgA, orange = sgB
#   so 3 green + 2 orange = five bands = that species has 5 copies, 3 with sgA
#   grey strip above a chromosome row = Drosera GENESPACE block, shade alternates at boundaries
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ggnewscale)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
W <- 8
SPORD <- c("regia","binata","paradoxa","scorpioides","capensis")
ABB <- c(regia="reg", binata="bin", paradoxa="par", scorpioides="sco", capensis="cap")
GREY <- c("#cccccc", "#8a8a8a")
BH <- 0.63; SH0 <- 0.22; SH1 <- 0.45     # barcode height, grey strip band

blk <- read_csv(file.path(GSD,"syntenicBlock_coordinates.csv"), show_col_types = FALSE)
v <- read_csv("tract_votes_blocks7.csv", show_col_types = FALSE) %>%
  mutate(species = sub("Drosera_","",sp),
         bid  = sub(".*:\\s*","", blk),
         dchr = sub("_hap1$|_collapsed$","", sub("^chr","c", lin_chr)),
         unit = paste(ABB[species], dchr))

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

## ---- slots: one per anchor gene, fixed across all rows ---------------------
anch <- v %>% distinct(pair, chrA, chrB, nep_chr, anchor, posA)
anch$dioblk <- mapply(dbf, anch$chrA, anch$posA, anch$nep_chr)
anch <- anch %>% arrange(pair, posA) %>% group_by(pair) %>%
  mutate(idx = row_number(), win = (idx - 1) %/% W + 1) %>%
  group_by(pair, win) %>% mutate(slot = row_number(), nslot = n()) %>% ungroup()
wsum <- anch %>% group_by(pair, chrA, chrB, nep_chr, win) %>%
  summarise(mid_Mb = median(posA)/1e6, nslot = n(),
            dioblk = names(sort(table(dioblk), decreasing=TRUE))[1], .groups="drop") %>%
  arrange(pair, win)
cat(sprintf("anchors %d | windows %d | votes %d\n", nrow(anch), nrow(wsum), nrow(v)))

V <- v %>% left_join(select(anch, pair, anchor, win, slot, nslot), by = c("pair","anchor"))

## rows: chromosomes, then per-species consensus, then ALL
ROW <- bind_rows(
  V %>% transmute(pair, anchor, win, slot, nslot, species, unit, vote, bid, kind="chr"),
  V %>% transmute(pair, anchor, win, slot, nslot, species,
                  unit = paste0("~ ", ABB[species], " tally"),
                  vote, bid = NA_character_, kind="spcons"),
  V %>% transmute(pair, anchor, win, slot, nslot, species = "consensus",
                  unit = "~~ ALL species", vote, bid = NA_character_, kind="allcons"),
  V %>% filter(species != "regia") %>%
    transmute(pair, anchor, win, slot, nslot, species = "consensus",
              unit = "~~ ALL minus regia", vote, bid = NA_character_, kind="allcons"))

## explicit display order: regia, binata, paradoxa, scorpioides, capensis,
## then the same order for the tallies, then ALL, then ALL minus regia
CONSORD <- c(paste0("~ ", ABB[SPORD], " tally"), "~~ ALL species", "~~ ALL minus regia")

pdf("FIG31_barcode_gene_aligned.pdf", width = 22, height = 12, onefile = TRUE)
np <- 0; SUMM <- list()
for (pp in sort(unique(ROW$pair))) {
  d <- filter(ROW, pair == pp); ws <- filter(wsum, pair == pp)
  if (!nrow(d)) next

  meta <- d %>% distinct(species, unit, kind) %>%
    left_join(d %>% group_by(unit) %>% summarise(fA = mean(vote=="A"), .groups="drop"),
              by="unit") %>%
    mutate(spf = factor(ifelse(kind=="chr", species, "zz_consensus"),
                        c(SPORD,"zz_consensus")),
           k1 = as.integer(spf),
           k2 = ifelse(kind == "chr", 0L, match(unit, CONSORD)),
           k3 = ifelse(kind == "chr", -fA, 0)) %>%
    filter(!is.na(spf)) %>% arrange(k1, k2, k3)
  d <- filter(d, unit %in% meta$unit)
  d$unit <- factor(d$unit, levels = rev(meta$unit))

  ## one band per gene copy, A first (bottom)
  dd <- d %>% arrange(unit, win, slot, vote) %>%
    group_by(unit, win, slot) %>% mutate(m = row_number(), k = n()) %>% ungroup() %>%
    mutate(yi = as.integer(unit), xi = match(win, ws$win)) %>% filter(!is.na(xi)) %>%
    mutate(sw = 0.92 / nslot,
           x0 = xi - 0.46 + (slot - 1) * sw + sw * 0.06,
           x1 = xi - 0.46 + slot * sw - sw * 0.06,
           bh = BH / k,
           y0 = yi - 0.45 + (m - 1) * bh,
           y1 = yi - 0.45 + m * bh - min(0.012, bh * 0.18),
           col = ifelse(vote == "A", "#1b9e77", "#d95f02"))

  ## grey block strip, chromosome rows only
  strip <- dd %>% filter(kind == "chr", !is.na(bid)) %>%
    distinct(unit, yi, win, slot, x0, x1, bid) %>%
    arrange(unit, win, slot) %>% group_by(unit) %>%
    mutate(chg = bid != lag(bid), chg = ifelse(is.na(chg), TRUE, chg),
           shade = GREY[(cumsum(chg) - 1) %% 2 + 1]) %>% ungroup()

  bb <- ws %>% mutate(chg = dioblk != lag(dioblk)) %>% filter(chg) %>% pull(win)
  rr <- meta %>% mutate(y = rev(seq_len(n())))
  div <- rr %>% group_by(spf) %>% summarise(top = max(y), .groups="drop") %>%
    mutate(cut = top + .5) %>% pull(cut)
  div <- div[div > .5 & div < nrow(rr)]
  lab_sp <- rr %>% group_by(spf) %>%
    summarise(y = mean(range(y)),
              lab = ifelse(spf[1]=="zz_consensus","CONSENSUS", as.character(spf[1])),
              .groups="drop")

  p <- ggplot() +
    geom_rect(data = strip, aes(xmin=x0, xmax=x1, ymin=yi+SH0, ymax=yi+SH1, fill=shade)) +
    scale_fill_identity(guide = "none") +
    new_scale_fill() +
    geom_rect(data = dd, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1, fill=vote)) +
    geom_rect(data = distinct(dd, unit, win, xi, yi),
              aes(xmin=xi-.47, xmax=xi+.47, ymin=yi-.47, ymax=yi+.47),
              fill=NA, colour="grey40", linewidth=.25) +
    { if (length(bb)) geom_vline(xintercept = match(bb, ws$win) - .5,
                                 linetype = 2, colour = "black", linewidth = .55) } +
    geom_hline(yintercept = div, colour = "black", linewidth = .8) +
    geom_text(data = lab_sp, aes(x = -1.3, y = y, label = lab),
              angle = 90, size = 3.6, fontface = 2, inherit.aes = FALSE) +
    scale_fill_manual(values = c(A="#1b9e77", B="#d95f02"), breaks = c("A","B"),
                      labels = c(A = sprintf("with sgA  (%s)", ws$chrA[1]),
                                 B = sprintf("with sgB  (%s)", ws$chrB[1])),
                      name = "gene tree places this Drosera copy") +
    scale_y_continuous(breaks = seq_len(nrow(rr)), labels = levels(d$unit),
                       expand = expansion(add = .7)) +
    scale_x_continuous(breaks = seq_along(ws$win),
                       labels = sprintf("%d\n%.0f", ws$win, ws$mid_Mb),
                       limits = c(-2.2, length(ws$win) + .6), expand = c(0,0)) +
    labs(title = sprintf("Nepenthes %s   |   Dionaea %s (sgA)  vs  %s (sgB)",
                         ws$nep_chr[1], ws$chrA[1], ws$chrB[1]),
         subtitle = paste0(
           "One SLOT = one anchor gene. Slots are fixed, so a slot can be read vertically across every species.\n",
           "One BAND inside a slot = one gene copy. Green bands stack at the bottom, orange on top: 3 green + 2 orange = 5 copies, 3 placed with sgA.\n",
           "Grey strip above a chromosome row = Drosera GENESPACE block; shade alternates at each boundary, so a colour flip inside one grey band means that block spans two ancestries.\n",
           sprintf("Columns = windows of %d consecutive anchor genes (window no. / Mb). Dashed lines = Dionaea block boundaries.", W)),
         x = "window along Dionaea sgA", y = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 6.5),
          axis.text.x = element_text(size = 6),
          legend.position = "top")
  print(p); np <- np + 1
  SUMM[[pp]] <- dd %>% filter(kind != "chr") %>%
    mutate(pair_lab = sprintf("%s  |  Dio %s vs %s", ws$nep_chr[1], ws$chrA[1], ws$chrB[1]))
}

## ---------- page 9: all consensus rows, every chromosome pair --------------
SU <- bind_rows(SUMM)
SU$unit <- factor(SU$unit, levels = rev(CONSORD))
SU <- SU %>% group_by(pair) %>% mutate(yi2 = as.integer(unit)) %>% ungroup() %>%
  arrange(pair, unit, win, slot, vote) %>%
  group_by(pair, unit, win, slot) %>% mutate(m = row_number(), k = n()) %>% ungroup() %>%
  mutate(bh = BH / k,
         y0 = yi2 - 0.45 + (m - 1) * bh,
         y1 = yi2 - 0.45 + m * bh - min(0.012, bh * 0.18))
ps <- ggplot(SU) +
  geom_rect(aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1, fill=vote)) +
  geom_rect(data = distinct(SU, pair_lab, unit, win, xi, yi2),
            aes(xmin=xi-.47, xmax=xi+.47, ymin=yi2-.47, ymax=yi2+.47),
            fill=NA, colour="grey40", linewidth=.2) +
  facet_wrap(~ pair_lab, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = c(A="#1b9e77", B="#d95f02"), breaks = c("A","B"),
                    labels = c(A="with sgA", B="with sgB"),
                    name = "gene tree places this Drosera copy") +
  scale_y_continuous(breaks = seq_along(CONSORD), labels = rev(CONSORD),
                     expand = expansion(add=.6)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(title = "SUMMARY: per-species tallies and pooled tallies, all Dionaea chromosome pairs",
       subtitle = paste0("One slot = one anchor gene, left to right along Dionaea sgA. One band = one gene copy.\n",
                         "Rows in fixed order: regia, binata, paradoxa, scorpioides, capensis, then ALL, then ALL minus regia.\n",
                         "sgA/sgB are LOCAL labels per pair and are not comparable between panels."),
       x = "window along Dionaea sgA", y = NULL) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(), legend.position = "top",
        axis.text.y = element_text(size = 6.5), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
print(ps); np <- np + 1
dev.off()

cat("\n=== slot occupancy: how often does a species have ANY copy of an anchor gene? ===\n")
occ <- V %>% distinct(pair, anchor, species) %>% count(pair, species, name="filled") %>%
  left_join(anch %>% count(pair, name="slots"), by="pair") %>%
  mutate(pct = round(100*filled/slots)) %>%
  select(pair, species, pct) %>% pivot_wider(names_from=species, values_from=pct)
print(as.data.frame(occ))
cat("(empty slots are real gene absences, not a plotting gap)\n")

cat("\n=== copies per gene per species (what the band stacks show) ===\n")
print(as.data.frame(V %>% count(pair, anchor, species, name="copies") %>%
  count(species, copies) %>% pivot_wider(names_from=copies, values_from=n, values_fill=0)))
cat(sprintf("\nWROTE: FIG31_barcode_gene_aligned.pdf (%d pages)\n", np))
