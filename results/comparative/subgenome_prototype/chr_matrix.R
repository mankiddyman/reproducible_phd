#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Rows = species x Drosera CHROMOSOME (no block-level filtering, nothing discarded).
# Cell  = stacked proportion of that chromosome's votes in that window (green sgA / orange sgB)
#         big number bottom-right = total votes
#         small text top-left     = GENESPACE block IDs contributing, with counts
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
W <- 8
SPORD <- c("regia","binata","paradoxa","scorpioides","capensis")
ABB <- c(regia="reg", binata="bin", paradoxa="par", scorpioides="sco", capensis="cap")

blk <- read_csv(file.path(GSD,"syntenicBlock_coordinates.csv"), show_col_types = FALSE)
v <- read_csv("tract_votes_blocks7.csv", show_col_types = FALSE) %>%
  mutate(species = sub("Drosera_","",sp),
         bid = sub(".*:\\s*","", blk),
         dchr = sub("_hap1$|_collapsed$","", sub("^chr","c", lin_chr)),
         unit = paste(ABB[species], dchr))
stopifnot(nrow(v) > 0)
cat(sprintf("votes in: %d | species x chromosome rows: %d\n", nrow(v), n_distinct(paste(v$pair, v$unit))))

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
  summarise(n_anch = n(), mid_Mb = median(posA)/1e6,
            dioblk = names(sort(table(dioblk), decreasing=TRUE))[1], .groups="drop") %>%
  arrange(pair, win)
stopifnot(nrow(wsum) > 0)

V <- v %>% left_join(select(anch, pair, anchor, win), by = c("pair","anchor"))

## cells: species x Drosera chromosome x window
CELL <- V %>% group_by(pair, species, unit, win) %>%
  summarise(nv = n(), fA = mean(vote == "A"),
            blocks = paste(sprintf("%s:%d", names(sort(table(bid), decreasing=TRUE)),
                                   sort(table(bid), decreasing=TRUE)), collapse=" "),
            nblk = n_distinct(bid), .groups="drop") %>% mutate(kind = "chr")
SPC <- V %>% group_by(pair, species, win) %>%
  summarise(nv = n(), fA = mean(vote == "A"), blocks = "", nblk = NA_integer_, .groups="drop") %>%
  mutate(unit = paste0("~ ", ABB[species], " consensus"), kind = "spcons")
ALLC <- SPC %>% group_by(pair, win) %>%
  summarise(nv = n(), fA = mean(fA), blocks = "", nblk = NA_integer_, .groups="drop") %>%
  mutate(species = "consensus", unit = "~~ ALL (mean of species)", kind = "allcons")
D <- bind_rows(CELL, SPC, ALLC)
stopifnot(nrow(D) > 0)

cat(sprintf("cells: %d chr-rows, %d species-consensus, %d ALL | votes displayed %d of %d (%.0f%%)\n",
            sum(D$kind=="chr"), sum(D$kind=="spcons"), sum(D$kind=="allcons"),
            sum(CELL$nv), nrow(v), 100*sum(CELL$nv)/nrow(v)))
cat("\nrows per species per Dionaea pair:\n")
print(as.data.frame(CELL %>% distinct(pair, species, unit) %>% count(pair, species) %>%
  pivot_wider(names_from=species, values_from=n, values_fill=0)))
cat("\nblock fragmentation per cell:\n"); print(summary(CELL$nblk))

pdf("FIG29_chr_matrix.pdf", width = 17, height = 11, onefile = TRUE)
np <- 0
for (pp in sort(unique(D$pair))) {
  d <- filter(D, pair == pp); ws <- filter(wsum, pair == pp)
  if (!sum(d$kind == "chr")) next

  rows <- d %>% distinct(species, unit, kind) %>%
    left_join(d %>% group_by(unit) %>% summarise(m = weighted.mean(fA, nv), .groups="drop"),
              by = "unit") %>%
    mutate(spf = factor(ifelse(kind == "chr", species, "zz"), c(SPORD, "zz"))) %>%
    filter(!is.na(spf)) %>% arrange(spf, kind, desc(m))
  d <- filter(d, unit %in% rows$unit)
  d$unit <- factor(d$unit, levels = rev(rows$unit))
  rr <- rows %>% mutate(y = rev(seq_len(n())))
  div <- rr %>% group_by(spf) %>% summarise(top = max(y), .groups="drop") %>% pull(top) + .5
  div <- div[div > .5 & div < nrow(rr)]
  lsp <- rr %>% group_by(spf) %>%
    summarise(y = mean(range(y)),
              lab = ifelse(spf[1]=="zz","CONSENSUS", as.character(spf[1])), .groups="drop")

  dd <- d %>% mutate(yi = as.integer(unit), xi = match(win, ws$win)) %>% filter(!is.na(xi)) %>%
    mutate(xmin = xi-.47, xmax = xi+.47,
           y0 = yi-.46, y1 = yi-.46+.92*fA, y2 = yi+.46)
  bb <- ws %>% mutate(chg = dioblk != lag(dioblk)) %>% filter(chg) %>% pull(win)

  p <- ggplot(dd) +
    geom_rect(aes(xmin=xmin, xmax=xmax, ymin=y1, ymax=y2), fill="#d95f02") +
    geom_rect(aes(xmin=xmin, xmax=xmax, ymin=y0, ymax=y1), fill="#1b9e77") +
    geom_rect(aes(xmin=xmin, xmax=xmax, ymin=y0, ymax=y2), fill=NA,
              colour="grey40", linewidth=.25) +
    geom_text(aes(xmin+.03, y2-.06, label = blocks), hjust=0, vjust=1,
              size=1.5, colour="grey20", lineheight=.8) +
    geom_text(aes(xmax-.04, y0+.06, label = nv), hjust=1, vjust=0,
              size=2.4, fontface=2, colour="grey5") +
    { if (length(bb)) geom_vline(xintercept = match(bb, ws$win)-.5,
                                 linetype=2, colour="black", linewidth=.5) } +
    geom_hline(yintercept = div, colour="black", linewidth=.7) +
    geom_text(data = lsp, aes(x=-1.3, y=y, label=lab), angle=90,
              size=3.4, fontface=2, inherit.aes=FALSE) +
    scale_y_continuous(breaks=seq_len(nrow(rr)), labels=levels(d$unit),
                       expand=expansion(add=.7)) +
    scale_x_continuous(breaks=seq_along(ws$win),
                       labels=sprintf("%d\n%.0f", ws$win, ws$mid_Mb),
                       limits=c(-2.2, length(ws$win)+.6), expand=c(0,0)) +
    labs(title = sprintf("Nepenthes %s   |   Dionaea %s (sgA) vs %s (sgB)",
                         ws$nep_chr[1], ws$chrA[1], ws$chrB[1]),
         subtitle = sprintf("row = species x Drosera chromosome, ALL votes used. Cell = stacked proportion: GREEN sgA (bottom) / ORANGE sgB (top).\nbold number bottom-right = votes; small text top-left = GENESPACE block IDs with counts; columns = %d-anchor windows (no. / Mb); dashed = Dionaea block boundary", W),
         x = "window along Dionaea sgA", y = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 7),
          axis.text.x = element_text(size = 6))
  print(p); np <- np + 1
}
dev.off()
write_csv(D, "chr_matrix_rows.csv")
cat(sprintf("\nWROTE: FIG29_chr_matrix.pdf (%d pages) chr_matrix_rows.csv\n", np))
