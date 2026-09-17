#!/usr/bin/env Rscript
# ============================================================================
# DR21 — COPY NUMBER PER LOCUS, NOT GENE CONTENT
#
# WHY DR19 AND DR20 BOTH FAILED TO TEST THE CLAIM
#   DR19 counted SEGMENTS. A segment is a contiguous run, so one ancestral copy
#   broken by a fusion counts twice. Segment counts measure fragmentation.
#
#   DR20 counted GENES. But A is DEFINED as the more-retained subgenome
#   (34-41% asymmetry in Dionaea), so
#       gene ratio = copy-number ratio  x  retention bias
#                  =        2.0         x     1.0 to 1.4
#   Gene content conflates the quantity being tested with a known bias. It
#   cannot test 2:1 at all. That was my design error, not a data problem.
#
# THIS TEST
#   Per locus, per species: how many tips are labelled A, how many B?
#   Fractionation only REMOVES copies, so each locus gives a FLOOR. Across
#   thousands of loci the maximum recovers true copy number, because at some
#   loci nothing was lost.
#
#     AAB hexaploid            -> max 2 A, 1 B
#     AAB doubled (12-ploid)   -> max 4 A, 2 B   (capensis, if DR08 is right)
#
#   No fractionation confound: we are counting copies, not surviving sequence.
#
# IN   DR/locus_meta.tsv, DR/out/DR02_propagated_blocks.csv
# OUT  DR/out/DR21_copynumber.csv
# FIG  DR/fig/DR21_1_copies.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

hr("1. LABEL EVERY TIP")
LOC <- read_tsv("DR/locus_meta.tsv", show_col_types=FALSE) %>%
       filter(genome %in% DROS)
PB <- read_csv("DR/out/DR02_propagated_blocks.csv", show_col_types=FALSE)
cat(sprintf("  Drosera tips in locus_meta: %d\n", nrow(LOC)))
cat(sprintf("  labelled genes available  : %d\n", nrow(PB)))
key <- setNames(PB$label, paste(PB$genome, PB$gene))
LOC$lab <- unname(key[paste(LOC$genome, LOC$gene)])
cat(sprintf("  tips with a label: %d (%.1f%%)\n",
            sum(!is.na(LOC$lab)), 100*mean(!is.na(LOC$lab))))
L <- LOC %>% filter(!is.na(lab))

hr("2. COPIES PER LOCUS")
CN <- L %>% count(genome, locus, lab) %>%
  pivot_wider(names_from=lab, values_from=n, values_fill=0) %>%
  mutate(sp=sub("Drosera_","",genome))
if (!"A" %in% names(CN)) CN$A <- 0L
if (!"B" %in% names(CN)) CN$B <- 0L
cat(sprintf("  loci x species with >=1 labelled copy: %d\n", nrow(CN)))

cat("\n  DISTRIBUTION of A copies per locus:\n")
print(CN %>% count(sp, A) %>% pivot_wider(names_from=A, values_from=n,
      values_fill=0, names_prefix="A=") %>% as.data.frame(), row.names=FALSE)
cat("\n  DISTRIBUTION of B copies per locus:\n")
print(CN %>% count(sp, B) %>% pivot_wider(names_from=B, values_from=n,
      values_fill=0, names_prefix="B=") %>% as.data.frame(), row.names=FALSE)

hr("3. THE CEILING — high quantiles recover true copy number")
cat("  Fractionation only removes copies, so the upper tail is the answer.\n\n")
q <- CN %>% group_by(sp) %>%
  summarise(loci=n(),
            A_p95=quantile(A,.95), A_p99=quantile(A,.99), A_max=max(A),
            B_p95=quantile(B,.95), B_p99=quantile(B,.99), B_max=max(B),
            .groups="drop")
print(as.data.frame(q), row.names=FALSE)
cat("\n  AAB predicts A_max 2, B_max 1.\n")
cat("  AAB doubled predicts A_max 4, B_max 2.\n")

hr("4. ONLY LOCI WHERE NOTHING WAS LOST")
cat("  Restrict to loci retaining the MOST copies -- these are the ones that\n")
cat("  escaped fractionation, so they show the underlying constitution.\n\n")
full <- CN %>% group_by(sp) %>%
  filter(A + B >= quantile(A+B, .90)) %>%
  summarise(n=n(), medA=median(A), medB=median(B),
            ratio=round(median(A)/pmax(median(B),1),2),
            modeAB=names(sort(table(paste0(A,"A:",B,"B")), decreasing=TRUE))[1],
            .groups="drop")
print(as.data.frame(full), row.names=FALSE)
cat("\n  modeAB is the single most common configuration among those loci.\n")

hr("5. VERDICT PER SPECIES")
for (g in unique(CN$sp)) {
  z <- CN[CN$sp==g, ]
  top <- sort(table(paste0(z$A,"A:",z$B,"B")), decreasing=TRUE)[1:3]
  cat(sprintf("\n  %s\n", g))
  cat(sprintf("    commonest configurations: %s\n",
              paste(names(top), top, sep=" x", collapse="  |  ")))
  cat(sprintf("    max copies seen: %dA / %dB\n", max(z$A), max(z$B)))
}
write_csv(CN, "DR/out/DR21_copynumber.csv")

PL <- CN %>% count(sp, A, B) %>% filter(A <= 8, B <= 6)
p1 <- ggplot(PL, aes(B, A, size=n, colour=sp)) +
  geom_abline(slope=2, intercept=0, linetype="dashed", colour="grey50") +
  geom_point(alpha=0.75) +
  facet_wrap(~sp) + scale_size_continuous(range=c(1,9), name="loci") +
  scale_x_continuous(breaks=0:6) + scale_y_continuous(breaks=0:8) +
  labs(title="DR21 - copies per locus, by subgenome",
       subtitle=paste0("dashed = 2:1. AAB predicts mass at (1,2); ",
                       "AAB doubled predicts mass at (2,4).\n",
                       "Fractionation pulls points toward the origin, ",
                       "so the OUTER envelope is what matters."),
       x="B copies at that locus", y="A copies at that locus", colour=NULL) +
  theme_minimal(10) + theme(legend.position="right")
suppressWarnings(ggsave("DR/fig/DR21_1_copies.pdf", p1, width=10, height=7))
cat("\n  wrote DR/fig/DR21_1_copies.pdf\n")
