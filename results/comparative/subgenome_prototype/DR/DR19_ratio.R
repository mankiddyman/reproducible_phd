#!/usr/bin/env Rscript
# ============================================================================
# DR19 — THE A:B RATIO WITHOUT PRE-SELECTION
#
# THE PROBLEM WITH THE 19/19 RESULT
#   It takes (species, region) pairs with EXACTLY 3 segments and asks whether
#   they are 2A:1B. But selecting on "has 3 segments" is close to selecting on
#   the answer: a hexaploid region IS three copies. 19 of 40 pairs qualify;
#   capensis contributes ZERO because it is 12-ploid and never has exactly 3.
#   So "single hybridisation in core Drosera" currently rests on four species,
#   none of them the 12-ploid.
#
# WHAT THIS DOES INSTEAD
#   Plots the A count against the B count for ALL 40 (species, region) pairs.
#   No filter. The claim predicts a line through the origin with slope 2 --
#   including capensis, which should sit near (4,2) rather than (2,1) but on
#   the SAME line. Ratio, not count, is the testable version of the claim.
#
# IN   DR/out/DR02_segments.csv
# OUT  DR/out/DR19_ratios.csv
# FIG  DR/fig/DR19_1_scatter.pdf, DR19_2_bars.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

S <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE) %>%
     filter(label %in% c("A","B"))
R <- S %>% count(genome, region, label) %>%
     pivot_wider(names_from=label, values_from=n, values_fill=0) %>%
     mutate(sp=sub("Drosera_","",genome), total=A+B,
            ratio=A/pmax(B,1), pattern=paste0(A,"A:",B,"B"))

hr("1. ALL 40 (species, region) PAIRS — nothing excluded")
print(as.data.frame(R %>% select(sp, region, A, B, total, ratio) %>%
      arrange(sp, region)), row.names=FALSE)

hr("2. THE OLD TEST vs AN UNFILTERED ONE")
cat(sprintf("  pairs with exactly 3 segments: %d of %d\n",
            sum(R$total==3), nrow(R)))
cat(sprintf("    of those, 2A:1B: %d\n", sum(R$total==3 & R$A==2 & R$B==1)))
cat("\n  UNFILTERED: is A > B in each pair?\n")
cat(sprintf("    A > B : %d\n    A = B : %d\n    A < B : %d\n",
            sum(R$A>R$B), sum(R$A==R$B), sum(R$A<R$B)))
bt <- binom.test(sum(R$A>R$B), sum(R$A!=R$B), 0.5)
cat(sprintf("    binomial p (A>B vs A<B, ties dropped) = %.3g\n", bt$p.value))
cat("\n  This uses all 40 pairs and does not condition on segment count.\n")

hr("3. IS THE RATIO 2, ACROSS ALL PAIRS?")
cat("  A 2:1 constitution predicts A = 2B. Fit through the origin.\n\n")
fit <- lm(A ~ 0 + B, data=R)
ci <- confint(fit)
cat(sprintf("  slope = %.3f  95%% CI [%.3f, %.3f]\n",
            coef(fit)[1], ci[1], ci[2]))
cat(sprintf("  does the CI contain 2? %s\n",
            ifelse(ci[1] <= 2 & ci[2] >= 2, "YES", "NO")))
cat("\n  per species:\n")
per <- R %>% group_by(sp) %>%
  summarise(pairs=n(), A=sum(A), B=sum(B), ratio=round(sum(A)/sum(B),2),
            .groups="drop")
print(as.data.frame(per), row.names=FALSE)
cat("\n  capensis is 12-ploid: expect higher COUNTS but the same RATIO.\n")

hr("4. WHICH PAIRS DEVIATE?")
R$expA <- 2*R$B
R$dev <- R$A - R$expA
cat("  pairs furthest from A = 2B:\n\n")
print(as.data.frame(R %>% arrange(desc(abs(dev))) %>%
      select(sp, region, A, B, expected_A=expA, deviation=dev) %>% head(10)),
      row.names=FALSE)
write_csv(R, "DR/out/DR19_ratios.csv")

p1 <- ggplot(R, aes(B, A, colour=sp)) +
  geom_abline(slope=2, intercept=0, linetype="dashed", colour="grey40") +
  geom_abline(slope=1, intercept=0, linetype="dotted", colour="grey70") +
  geom_point(size=3, alpha=0.85) +
  ggrepel::geom_text_repel(aes(label=sub("_dom","",region)), size=2.6,
                           show.legend=FALSE, max.overlaps=20) +
  facet_wrap(~sp, scales="free") +
  labs(title="DR19.1 - A count vs B count, all 40 (species, region) pairs",
       subtitle=paste0("dashed = the 2:1 prediction, dotted = 1:1. ",
         "No pair is excluded.\ncapensis is 12-ploid so its points sit ",
         "further out, but should lie on the SAME dashed line."),
       x="B segments", y="A segments", colour=NULL) +
  theme_minimal(9) + theme(legend.position="none")
suppressWarnings(ggsave("DR/fig/DR19_1_scatter.pdf", p1, width=10, height=7))

p2 <- ggplot(R, aes(reorder(paste(sp, sub("_dom","",region)), ratio), ratio,
                    fill=sp)) +
  geom_hline(yintercept=2, linetype="dashed", colour="grey30") +
  geom_col(width=0.7) + coord_flip() +
  labs(title="DR19.2 - A:B ratio per (species, region)",
       subtitle="dashed = 2:1. All 40 pairs, unfiltered.",
       x=NULL, y="A / B", fill=NULL) +
  theme_minimal(8) + theme(legend.position="top")
suppressWarnings(ggsave("DR/fig/DR19_2_bars.pdf", p2, width=8, height=10))
cat("\n  wrote DR/fig/DR19_1_scatter.pdf and DR19_2_bars.pdf\n")
