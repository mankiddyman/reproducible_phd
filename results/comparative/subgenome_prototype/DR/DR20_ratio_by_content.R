#!/usr/bin/env Rscript
# ============================================================================
# DR20 — A:B BY SEQUENCE CONTENT, NOT SEGMENT COUNT
#
# WHY DR19 SCATTERS
#   A segment is a contiguous run of one ancestral region on one chromosome.
#   A single ancestral COPY broken by a fusion appears as TWO segments. So
#   segment counts measure fragmentation as much as copy number, and the
#   species differ enormously in fragmentation: capensis 20 chromosomes,
#   binata 16, regia 17, paradoxa 6, scorpioides 4.
#
#   Counting Mb or genes is fragmentation-invariant: however a copy is broken
#   up, its total sequence is the same.
#
# THE THREE THINGS THIS SEPARATES
#   1. ratio by CONTENT vs by COUNT -- does the scatter shrink?
#   2. does deviation from 2 track VOTING GENE DENSITY? if the noisiest
#      regions deviate most, the scatter is detection noise, not biology
#   3. capensis on its own terms -- 12-ploid, so 4A:2B is still ratio 2
#
# IN   DR/out/DR02_propagated_blocks.csv, DR02_segments.csv, DR02_gene_delta.csv
# OUT  DR/out/DR20_content.csv
# FIG  DR/fig/DR20_1_content.pdf, DR20_2_noise.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

PB <- read_csv("DR/out/DR02_propagated_blocks.csv", show_col_types=FALSE)
S  <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE) %>%
      filter(label %in% c("A","B"))
GD <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE)

hr("1. RATIO BY GENE CONTENT vs BY SEGMENT COUNT")
byGene <- PB %>% count(genome, region, label) %>%
  pivot_wider(names_from=label, values_from=n, values_fill=0) %>%
  transmute(genome, region, gA=A, gB=B, gene_ratio=A/pmax(B,1))
bySeg <- S %>% count(genome, region, label) %>%
  pivot_wider(names_from=label, values_from=n, values_fill=0) %>%
  transmute(genome, region, sA=A, sB=B, seg_ratio=A/pmax(B,1))
R <- full_join(byGene, bySeg, by=c("genome","region")) %>%
  mutate(sp=sub("Drosera_","",genome))
print(as.data.frame(R %>% select(sp, region, sA, sB, seg_ratio, gA, gB, gene_ratio) %>%
      mutate(across(c(seg_ratio,gene_ratio), ~round(.,2))) %>% arrange(sp, region)),
      row.names=FALSE)

hr("2. WHICH MEASURE IS TIGHTER AROUND 2?")
f <- function(v) sprintf("median %.2f | IQR %.2f-%.2f | MAD from 2 = %.2f",
                         median(v), quantile(v,.25), quantile(v,.75),
                         median(abs(v-2)))
cat("  by SEGMENT COUNT :", f(R$seg_ratio), "\n")
cat("  by GENE CONTENT  :", f(R$gene_ratio), "\n")
cat("\n  If gene content is tighter, the scatter in DR19 was fragmentation.\n")

hr("3. IS DEVIATION JUST NOISE?")
cat("  If the regions with fewest voting genes deviate most, the spread is\n")
cat("  detection noise rather than real variation in constitution.\n\n")
dens <- GD %>% count(genome, region, name="nvote")
R <- R %>% left_join(dens, by=c("genome","region")) %>%
     mutate(dev = abs(gene_ratio - 2))
ct <- cor.test(R$nvote, R$dev, method="spearman")
cat(sprintf("  Spearman(voting genes, |ratio - 2|) = %.3f, p = %.3g\n",
            ct$estimate, ct$p.value))
cat(sprintf("  regions with <200 voting genes: median |dev| = %.2f\n",
            median(R$dev[R$nvote < 200])))
cat(sprintf("  regions with >=200            : median |dev| = %.2f\n",
            median(R$dev[R$nvote >= 200])))

hr("4. PER SPECIES, POOLED ACROSS REGIONS")
per <- R %>% group_by(sp) %>%
  summarise(regions=n(), genes_A=sum(gA), genes_B=sum(gB),
            pooled_ratio=round(sum(gA)/sum(gB),2),
            median_region_ratio=round(median(gene_ratio),2),
            .groups="drop")
print(as.data.frame(per), row.names=FALSE)
cat("\n  AAB predicts 2.0 for all five, including capensis: a 12-ploid from\n")
cat("  doubling an AAB hexaploid is 4A:2B, which is still ratio 2.\n")

hr("5. DOES THE POOLED RATIO DIFFER FROM 2?")
for (g in unique(R$sp)) {
  z <- R[R$sp==g, ]
  bt <- binom.test(sum(z$gA), sum(z$gA)+sum(z$gB), 2/3)
  cat(sprintf("  %-12s A=%6d B=%6d  ratio %.2f  p(vs 2:1) = %.3g\n",
              g, sum(z$gA), sum(z$gB), sum(z$gA)/sum(z$gB), bt$p.value))
}
cat("\n  NOTE: genes are not independent, so these p-values are\n")
cat("  anticonservative. Use them to rank species, not as evidence.\n")
write_csv(R, "DR/out/DR20_content.csv")

L <- R %>% select(sp, region, seg_ratio, gene_ratio) %>%
     pivot_longer(c(seg_ratio, gene_ratio), names_to="measure", values_to="ratio")
p1 <- ggplot(L, aes(sp, ratio, colour=measure)) +
  geom_hline(yintercept=2, linetype="dashed", colour="grey30") +
  geom_point(position=position_jitterdodge(jitter.width=0.15, dodge.width=0.6),
             size=2.4, alpha=0.8) +
  scale_colour_manual(values=c(seg_ratio="#D85A30", gene_ratio="#1D9E75"),
                      labels=c("gene content","segment count"), name=NULL) +
  labs(title="DR20.1 - A:B by gene content vs by segment count",
       subtitle="all 40 (species, region) pairs, nothing excluded | dashed = 2:1",
       x=NULL, y="A / B") +
  theme_minimal(10) + theme(legend.position="top")
suppressWarnings(ggsave("DR/fig/DR20_1_content.pdf", p1, width=9, height=6))

p2 <- ggplot(R, aes(nvote, gene_ratio, colour=sp)) +
  geom_hline(yintercept=2, linetype="dashed", colour="grey30") +
  geom_point(size=2.6) + scale_x_log10() +
  labs(title="DR20.2 - does deviation from 2:1 track data volume?",
       subtitle="if the low-coverage regions scatter most, the spread is noise",
       x="voting genes in that region (log)", y="A / B by gene content",
       colour=NULL) +
  theme_minimal(10) + theme(legend.position="top")
suppressWarnings(ggsave("DR/fig/DR20_2_noise.pdf", p2, width=9, height=6))
cat("\n  wrote DR/fig/DR20_1_content.pdf and DR20_2_noise.pdf\n")
