#!/usr/bin/env Rscript
# ============================================================================
# DR15 — CHROMOSOME PAINTING: WHERE ARE A AND B?
#
# The GENESPACE riparian shows SYNTENY -- ribbons between chromosomes, coloured
# by the reference genome's chromosome. It cannot show subgenome assignment,
# because that is not in the synteny data.
#
# This uses the same colour scheme (ancestral region = Nepenthes _dom
# chromosome, red through blue) but paints the assignment onto each Drosera
# chromosome directly:
#
#   each chromosome is a bar
#   A segments sit in the UPPER half, B segments in the LOWER half
#   fill colour = which ancestral region the segment came from
#   black outline = a MARGINAL segment (vote purity < 0.6)
#
# So an AAB hexaploid should show, for each colour, TWO chromosomes painted on
# the A track and ONE on the B track. The constitution becomes visible rather
# than tabular.
#
# WHAT THE FIGURE HIDES, AND WHY IT IS ACCEPTABLE
#   Only 8-19% of genes have a directly measured delta; the rest inherit their
#   segment's label. DR14 showed 94.6-98.3% of propagated genes sit within
#   100 kb of a measured one, so this is dense interpolation rather than
#   extrapolation. The 0.1-1.2% beyond 5 Mb are the exception.
#
#   Vote purity averages 0.68-0.74. That is close to the CEILING: per-gene
#   delta has ~0.5 SD separation, capping per-gene accuracy near 0.60. What
#   makes the calls reliable is vote COUNT -- at 0.7 purity and 50 votes,
#   P(correct) = 0.998, and the median segment has 44-71 votes.
#
#   The 12 segments below 0.6 purity are outlined in black. Seven are regia.
#
# IN   DR/out/DR02_segments.csv, DR14_segment_quality.csv,
#      ../genespace/results/combBed.txt
# OUT  DR/fig/DR15_1_painting.pdf, DR15_2_constitution.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))
GSD <- file.path(dirname(getwd()), "genespace", "results")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

hr("1. LOAD")
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE) %>%
      filter(label %in% c("A","B"))
PUR <- read_csv("DR/out/DR14_segment_quality.csv", show_col_types=FALSE) %>%
       select(segment, purity, votes)
SS <- SS %>% left_join(PUR, by="segment") %>%
      mutate(marginal = !is.na(purity) & purity < 0.6)
cat(sprintf("  segments: %d | marginal: %d\n", nrow(SS), sum(SS$marginal)))

bed <- read.table(file.path(GSD,"combBed.txt"), header=TRUE, sep="\t",
                  quote="", comment.char="", stringsAsFactors=FALSE)
CHRLEN <- bed %>% group_by(genome, chr) %>%
  summarise(len_mb = max(end)/1e6, genes = n(), .groups="drop") %>%
  filter(genes >= 50)
cat(sprintf("  chromosomes with >=50 genes: %d\n", nrow(CHRLEN)))

# match the riparian palette: red -> blue across the Nepenthes _dom regions
REG <- sort(unique(SS$region))
cat("\n  ancestral regions:", paste(REG, collapse=", "), "\n")
PAL <- setNames(colorRampPalette(
  c("#C0392B","#E67E22","#F1C40F","#D4E157","#AED581",
    "#4FC3F7","#2980B9","#1A5276"))(length(REG)), REG)

hr("2. WHAT GETS DRAWN")
plotd <- SS %>%
  inner_join(CHRLEN, by=c("genome","chr")) %>%
  mutate(sp = sub("Drosera_","",genome),
         ymin = ifelse(label=="A", 0.05, -0.45),
         ymax = ifelse(label=="A", 0.45, -0.05))
print(as.data.frame(plotd %>% count(sp, label) %>%
  pivot_wider(names_from=label, values_from=n, values_fill=0)), row.names=FALSE)
cat("\n  chromosomes carrying at least one segment, per species:\n")
print(as.data.frame(plotd %>% group_by(sp) %>%
  summarise(chromosomes=n_distinct(chr), segments=n(), .groups="drop")),
  row.names=FALSE)

hr("3. THE PAINTING")
bg <- plotd %>% distinct(sp, chr, len_mb)
p <- ggplot() +
  geom_rect(data=bg, aes(xmin=0, xmax=len_mb, ymin=-0.45, ymax=0.45),
            fill="grey92", colour="grey75", linewidth=0.2) +
  geom_hline(yintercept=0, colour="white", linewidth=0.4) +
  geom_rect(data=plotd,
            aes(xmin=mb_lo, xmax=mb_hi, ymin=ymin, ymax=ymax, fill=region),
            colour=NA) +
  geom_rect(data=plotd %>% filter(marginal),
            aes(xmin=mb_lo, xmax=mb_hi, ymin=ymin, ymax=ymax),
            fill=NA, colour="black", linewidth=0.45) +
  scale_fill_manual(values=PAL, name="ancestral region\n(Nepenthes chromosome)") +
  facet_grid(chr ~ sp, scales="free", space="free_y", switch="y") +
  labs(title="Subgenome assignment along Drosera chromosomes",
       subtitle=paste0("UPPER half of each bar = A subgenome, LOWER half = B. ",
         "Colour = ancestral region, matching the riparian.\n",
         "Black outline = marginal segment (vote purity < 0.6, 12 of 190)."),
       x="position (Mb)", y=NULL) +
  theme_minimal(8) +
  theme(strip.text.y.left = element_text(angle=0, hjust=1, size=5),
        strip.text.x = element_text(face="bold", size=9),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "bottom",
        plot.subtitle = element_text(size=8, colour="grey30"))
ggsave("DR/fig/DR15_1_painting.pdf", p, width=13, height=16, limitsize=FALSE)
cat("  wrote DR/fig/DR15_1_painting.pdf\n")

hr("4. CONSTITUTION, drawn")
cat("  For an AAB hexaploid: two A segments and one B per ancestral region.\n\n")
CON <- SS %>% count(genome, region, label) %>%
  mutate(sp = sub("Drosera_","",genome)) %>%
  complete(nesting(genome, sp, region), label=c("A","B"), fill=list(n=0L))
print(as.data.frame(CON %>% select(sp, region, label, n) %>%
  pivot_wider(names_from=label, values_from=n, values_fill=0) %>%
  mutate(pattern=paste0(A,"A:",B,"B")) %>%
  arrange(sp, region)), row.names=FALSE)
p2 <- ggplot(CON, aes(region, n, fill=label)) +
  geom_col(position="dodge", width=0.7) +
  geom_hline(yintercept=c(1,2), linetype="dotted", colour="grey40") +
  facet_wrap(~sp, ncol=2) +
  scale_fill_manual(values=c(A="#1D9E75", B="#D85A30")) +
  labs(title="Segment counts per ancestral region",
       subtitle="AAB predicts 2 A and 1 B per region | dotted lines at 1 and 2",
       x="ancestral region", y="segments", fill="subgenome") +
  theme_minimal(9) +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top")
ggsave("DR/fig/DR15_2_constitution.pdf", p2, width=10, height=8)
cat("  wrote DR/fig/DR15_2_constitution.pdf\n")

hr("DONE")
cat("  FIG 1 -- where A and B sit along each chromosome\n")
cat("  FIG 2 -- whether the 2A:1B constitution holds region by region\n")
