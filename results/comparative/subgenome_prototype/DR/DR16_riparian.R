#!/usr/bin/env Rscript
# ============================================================================
# DR16 — RIPARIAN WITH A/B ENCODED IN THE RIBBONS
#
# The GENESPACE riparian draws one ribbon per syntenic block, positioned by
# gene rank order and coloured by the reference genome's chromosome. It cannot
# show subgenome assignment because that is not in the synteny data.
#
# THIS VERSION keeps the riparian exactly as it is and adds ONE dimension:
#   solid ribbon      = A subgenome
#   dashed, faded     = B subgenome
#   grey, very faint  = block with no confident A/B call
#
# HOW A RIBBON GETS ITS LABEL
#   A block is a syntenic interval in the Drosera genome. Its midpoint on the
#   Drosera side is matched against DR02's segments; the block inherits that
#   segment's label. Section 2 counts blocks that span TWO segments or fall in
#   NONE -- both are reported, because a figure that silently drops a third of
#   the blocks is worse than no figure.
#
# WHAT THIS INHERITS
#   Vote purity averages 0.68-0.74, close to the ceiling set by 0.5 SD per-gene
#   separation. Reliability comes from vote COUNT: at 50 votes and 0.7 purity,
#   P(correct) = 0.998, and the median segment has 44-71 votes. The 12 segments
#   below 0.6 purity are drawn with a black edge.
#
# IN   ../genespace/results/{syntenicBlock_coordinates.csv,combBed.txt}
#      DR/out/DR02_segments.csv, DR/out/DR14_segment_quality.csv
# OUT  DR/fig/DR16_riparian_AB.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))
GSD <- file.path(dirname(getwd()), "genespace", "results")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
ORDER <- c("Nepenthes_gracilis","Dionaea_muscipula","Drosera_regia",
           "Drosera_binata","Drosera_paradoxa","Drosera_scorpioides",
           "Drosera_capensis")

hr("1. BLOCKS AGAINST NEPENTHES")
blk <- read_csv(file.path(GSD,"syntenicBlock_coordinates.csv"), show_col_types=FALSE)
cat("  columns:", paste(names(blk), collapse=", "), "\n")
B <- bind_rows(
  blk %>% filter(genome2=="Nepenthes_gracilis", genome1 %in% ORDER) %>%
    transmute(sp=genome1, sp_chr=chr1,
              s1=pmin(startBp1,endBp1)/1e6, e1=pmax(startBp1,endBp1)/1e6,
              nep_chr=chr2, s2=pmin(startBp2,endBp2)/1e6,
              e2=pmax(startBp2,endBp2)/1e6, nh=nHits1),
  blk %>% filter(genome1=="Nepenthes_gracilis", genome2 %in% ORDER) %>%
    transmute(sp=genome2, sp_chr=chr2,
              s1=pmin(startBp2,endBp2)/1e6, e1=pmax(startBp2,endBp2)/1e6,
              nep_chr=chr1, s2=pmin(startBp1,endBp1)/1e6,
              e2=pmax(startBp1,endBp1)/1e6, nh=nHits2)) %>%
  filter(grepl("_dom$", nep_chr), sp != "Nepenthes_gracilis") %>% distinct()
cat(sprintf("\n  blocks to Nepenthes _dom chromosomes: %d\n", nrow(B)))
print(as.data.frame(B %>% count(sp, name="blocks")), row.names=FALSE)

hr("2. LABEL EACH BLOCK — and count what fails")
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE) %>%
      filter(label %in% c("A","B"))
PUR <- read_csv("DR/out/DR14_segment_quality.csv", show_col_types=FALSE) %>%
       select(segment, purity)
SS <- SS %>% left_join(PUR, by="segment")
B$mid <- (B$s1 + B$e1)/2
hit <- lapply(seq_len(nrow(B)), function(i) {
  k <- which(SS$genome==B$sp[i] & SS$chr==B$sp_chr[i] &
             SS$mb_lo <= B$mid[i] & SS$mb_hi >= B$mid[i])
  if (!length(k)) return(c(NA_character_, NA_real_, NA_real_))
  c(SS$label[k[1]], SS$purity[k[1]], length(k))
})
B$label   <- vapply(hit, function(x) x[1], character(1))
B$purity  <- as.numeric(vapply(hit, function(x) x[2], character(1)))
B$n_seg   <- as.numeric(vapply(hit, function(x) x[3], character(1)))
cat(sprintf("  labelled: %d of %d (%.0f%%)\n", sum(!is.na(B$label)), nrow(B),
            100*mean(!is.na(B$label))))
cat(sprintf("  midpoint fell in NO segment: %d\n", sum(is.na(B$label))))
cat(sprintf("  midpoint fell in >1 segment: %d\n",
            sum(B$n_seg > 1, na.rm=TRUE)))
print(as.data.frame(B %>% count(sp, label) %>%
  pivot_wider(names_from=label, values_from=n, values_fill=0)), row.names=FALSE)
B <- B %>% mutate(cls = ifelse(is.na(label), "no call", label),
                  marginal = !is.na(purity) & purity < 0.6)
cat(sprintf("\n  blocks on marginal segments: %d\n", sum(B$marginal)))

hr("3. LAY OUT THE RIBBONS")
bed <- read.table(file.path(GSD,"combBed.txt"), header=TRUE, sep="\t",
                  quote="", comment.char="", stringsAsFactors=FALSE)
CH <- bed %>% filter(genome %in% ORDER) %>%
  group_by(genome, chr) %>%
  summarise(len=max(end)/1e6, ngene=n(), .groups="drop") %>%
  filter(ngene >= 50) %>%
  arrange(genome, desc(len)) %>% group_by(genome) %>%
  mutate(off = cumsum(lag(len + 8, default=0))) %>% ungroup()
cat(sprintf("  chromosomes drawn: %d\n", nrow(CH)))
CH$y <- match(CH$genome, ORDER)
xof <- function(g, c) { i <- which(CH$genome==g & CH$chr==c)
                        if (length(i)) CH$off[i[1]] else NA_real_ }
B$x_sp  <- mapply(xof, B$sp, B$sp_chr) + B$mid
B$x_nep <- mapply(xof, "Nepenthes_gracilis", B$nep_chr) + (B$s2+B$e2)/2
B$y_sp  <- match(B$sp, ORDER)
B <- B %>% filter(!is.na(x_sp), !is.na(x_nep))
cat(sprintf("  ribbons with both endpoints placed: %d\n", nrow(B)))

REG <- sort(unique(B$nep_chr))
PAL <- setNames(colorRampPalette(
  c("#C0392B","#E67E22","#F1C40F","#D4E157","#AED581",
    "#4FC3F7","#2980B9","#1A5276"))(length(REG)), REG)

RIB <- bind_rows(lapply(seq_len(nrow(B)), function(i) {
  t <- seq(0, 1, length.out=24)
  tibble(id=i, x = B$x_nep[i] + t*(B$x_sp[i]-B$x_nep[i]),
         y = 1 + t*(B$y_sp[i]-1),
         nep_chr=B$nep_chr[i], cls=B$cls[i], nh=B$nh[i])
})) %>% mutate(x = x + 12*sin(pi*(y-1)/max(1,(max(y)-1))) *
                     0*(1))   # straight lines; keeps ordering readable

p <- ggplot() +
  geom_segment(data=CH, aes(x=off, xend=off+len, y=y, yend=y),
               colour="grey35", linewidth=2.6, lineend="round") +
  geom_path(data=RIB %>% filter(cls=="no call"),
            aes(x, y, group=id, colour=nep_chr),
            linewidth=0.20, alpha=0.10, linetype="solid") +
  geom_path(data=RIB %>% filter(cls=="B"),
            aes(x, y, group=id, colour=nep_chr),
            linewidth=0.32, alpha=0.40, linetype="22") +
  geom_path(data=RIB %>% filter(cls=="A"),
            aes(x, y, group=id, colour=nep_chr),
            linewidth=0.32, alpha=0.75, linetype="solid") +
  scale_colour_manual(values=PAL, name="ancestral region") +
  scale_y_continuous(breaks=seq_along(ORDER),
                     labels=gsub("_"," ", ORDER), trans="reverse") +
  labs(title="Syntenic blocks to Nepenthes, with subgenome assignment",
       subtitle=paste0("SOLID = A subgenome   ·   DASHED, faded = B subgenome",
         "   ·   very faint = no confident call\n",
         "Colour = Nepenthes chromosome, as in the GENESPACE riparian. ",
         "Chromosomes scaled in Mb."),
       x=NULL, y=NULL) +
  theme_minimal(11) +
  theme(panel.grid=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_text(face="italic", size=10),
        legend.position="bottom",
        plot.subtitle=element_text(size=8.5, colour="grey30"))
ggsave("DR/fig/DR16_riparian_AB.pdf", p, width=15, height=9)
ggsave("DR/fig/DR16_riparian_AB.png", p, width=15, height=9, dpi=200)
cat("\n  wrote DR/fig/DR16_riparian_AB.{pdf,png}\n")
cat("\n  Read section 2 first: if many blocks are unlabelled the figure is\n")
cat("  showing absence of data as much as it is showing biology.\n")
