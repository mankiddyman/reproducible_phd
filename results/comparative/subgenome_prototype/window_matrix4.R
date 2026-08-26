#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Multi-page PDF, one page per Dionaea chromosome pair.
#   columns = position-ordered windows of consecutive anchor genes
#   rows    = Drosera lineage blocks, grouped by species, horizontal dividers between groups
#   cell    = stacked proportion of that block's votes in that window (green sgA / orange sgB)
#   species consensus = pooled over EXACTLY the displayed blocks of that species
#   ALL consensus     = unweighted mean of the species fractions (copy number cannot dominate)
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
         bid  = sub(".*:\\s*", "", blk),
         unit = sprintf("%s %s b%s", ABB[species],
                        sub("_hap1$|_collapsed$","", sub("^chr","c",lin_chr)), bid))

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

## ---- anchors -> windows, with a uniqueness check --------------------------
anch <- v %>% distinct(pair, chrA, chrB, nep_chr, anchor, posA)
dup <- anch %>% count(pair, anchor) %>% filter(n > 1)
if (nrow(dup)) { cat("WARNING: anchors with >1 row (would corrupt windowing):\n"); print(as.data.frame(head(dup))) }
anch$dioblk <- mapply(dbf, anch$chrA, anch$posA, anch$nep_chr)
anch <- anch %>% arrange(pair, posA) %>% group_by(pair) %>%
  mutate(win = (row_number()-1) %/% W + 1) %>% ungroup()
wsum <- anch %>% group_by(pair, chrA, chrB, nep_chr, win) %>%
  summarise(n_anch = n(), mid_Mb = median(posA)/1e6,
            dioblk = names(sort(table(dioblk), decreasing=TRUE))[1], .groups="drop") %>%
  arrange(pair, win)

V <- v %>% left_join(select(anch, pair, anchor, win), by = c("pair","anchor"))

## ---- 1. decide which blocks are DISPLAYED ---------------------------------
keep <- V %>% count(pair, species, unit, win, name="k") %>%
  group_by(pair, unit) %>%
  filter(sum(k) >= MINROW, n_distinct(win) >= 2) %>% ungroup() %>%
  distinct(pair, species, unit)
Vk <- V %>% semi_join(keep, by = c("pair","species","unit"))
cat(sprintf("votes total %d | in displayed blocks %d (%.0f%%)\n",
            nrow(V), nrow(Vk), 100*nrow(Vk)/nrow(V)))

## ---- 2. rows, consensus derived from the SAME set --------------------------
BLK <- Vk %>% group_by(pair, win, species, unit) %>%
  summarise(nv = n(), fA = mean(vote == "A"), .groups="drop") %>% mutate(kind = "block")
SPC <- Vk %>% group_by(pair, win, species) %>%
  summarise(nv = n(), fA = mean(vote == "A"), .groups="drop") %>%
  mutate(unit = paste0("~ ", ABB[species], " consensus"), kind = "spcons")
ALLC <- SPC %>% group_by(pair, win) %>%
  summarise(nv = n(), fA = mean(fA), .groups="drop") %>%          # unweighted mean of species
  mutate(species = "consensus", unit = "~~ ALL (mean of species)", kind = "allcons")
D <- bind_rows(BLK, SPC, ALLC)

## ---- 3. reconciliation check ----------------------------------------------
chk <- BLK %>% group_by(pair, win, species) %>%
  summarise(sum_blocks = sum(nv), .groups="drop") %>%
  left_join(select(SPC, pair, win, species, cons_nv = nv), by=c("pair","win","species")) %>%
  mutate(ok = sum_blocks == cons_nv)
cat(sprintf("reconciliation: %d/%d (pair,window,species) cells where displayed blocks sum to the species consensus\n",
            sum(chk$ok), nrow(chk)))
if (any(!chk$ok)) print(as.data.frame(filter(chk, !ok) %>% head(10)))

cat("\ndisplayed blocks per species per Dionaea pair:\n")
print(as.data.frame(keep %>% count(pair, species) %>%
        pivot_wider(names_from=species, values_from=n, values_fill=0)))

## ---- 4. draw ---------------------------------------------------------------
pdf("FIG26_window_matrix_all.pdf", width = 16, height = 11, onefile = TRUE)
np <- 0
for (pp in sort(unique(D$pair))) {
  d <- filter(D, pair == pp); ws <- filter(wsum, pair == pp)
  if (sum(d$kind == "block") == 0) { cat(sprintf("skipping %s: no displayed blocks\n", pp)); next }

  rows <- d %>% distinct(species, unit, kind) %>%
    left_join(d %>% group_by(unit) %>% summarise(m = weighted.mean(fA, nv), .groups="drop"),
              by="unit") %>%
    mutate(spf = factor(ifelse(kind == "block", species, "zz_consensus"),
                        c(SPORD, "zz_consensus"))) %>%
    filter(!is.na(spf)) %>% arrange(spf, kind, desc(m))
  d <- filter(d, unit %in% rows$unit)
  d$unit <- factor(d$unit, levels = rev(rows$unit))

  rr <- rows %>% mutate(y = rev(seq_len(n())))
  div <- rr %>% group_by(spf) %>% summarise(top = max(y), .groups="drop") %>%
    mutate(cut = top + 0.5) %>% pull(cut)
  div <- div[div > 0.5 & div < nrow(rr)]
  labs_sp <- rr %>% group_by(spf) %>%
    summarise(y = mean(range(y)),
              lab = ifelse(spf[1] == "zz_consensus", "CONSENSUS", as.character(spf[1])),
              .groups="drop")

  dd <- d %>% mutate(yi = as.integer(unit), xi = match(win, ws$win)) %>%
    filter(!is.na(xi)) %>%
    mutate(xmin = xi - .46, xmax = xi + .46,
           yA0 = yi - .45, yA1 = yi - .45 + .9*fA, yB1 = yi + .45)

  bb <- ws %>% mutate(chg = dioblk != lag(dioblk)) %>% filter(chg) %>% pull(win)
  xl <- sprintf("%d\n%.0f", ws$win, ws$mid_Mb)

  p <- ggplot(dd) +
    geom_rect(aes(xmin=xmin, xmax=xmax, ymin=yA1, ymax=yB1), fill = "#d95f02") +
    geom_rect(aes(xmin=xmin, xmax=xmax, ymin=yA0, ymax=yA1), fill = "#1b9e77") +
    geom_rect(aes(xmin=xmin, xmax=xmax, ymin=yA0, ymax=yB1),
              fill=NA, colour="grey45", linewidth=.25) +
    geom_text(aes(xi, yi, label = nv), size = 2, colour = "grey10") +
    { if (length(bb)) geom_vline(xintercept = match(bb, ws$win) - .5,
                                 linetype = 2, colour = "black", linewidth = .5) } +
    geom_hline(yintercept = div, colour = "black", linewidth = .7) +
    geom_text(data = labs_sp, aes(x = -1.4, y = y, label = lab),
              angle = 90, size = 3.4, fontface = 2, inherit.aes = FALSE) +
    scale_y_continuous(breaks = seq_len(nrow(rr)), labels = levels(d$unit),
                       expand = expansion(add = .7)) +
    scale_x_continuous(breaks = seq_along(ws$win), labels = xl,
                       limits = c(-2.4, length(ws$win) + .6), expand = c(0,0)) +
    labs(title = sprintf("Nepenthes %s   |   Dionaea %s (sgA) vs %s (sgB)",
                         ws$nep_chr[1], ws$chrA[1], ws$chrB[1]),
         subtitle = sprintf("cell = stacked proportion of that block's votes in the window: GREEN = sgA (bottom), ORANGE = sgB (top); number = votes\ncolumns = %d-anchor windows left to right (window no. / Mb); dashed = GENESPACE Dionaea block boundary; ALL row = unweighted mean of the species consensus rows", W),
         x = "window along Dionaea sgA", y = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 6.5),
          axis.text.x = element_text(size = 6))
  print(p); np <- np + 1
}
dev.off()
write_csv(D, "window_matrix_rows.csv"); write_csv(wsum, "window_index2.csv")
cat(sprintf("\nWROTE: FIG26_window_matrix_all.pdf (%d pages)\n", np))
