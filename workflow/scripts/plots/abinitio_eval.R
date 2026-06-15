#!/usr/bin/env Rscript
# Ab initio evaluation (binata vs paradoxa), braker_filtered as RNA reference.
# Fig A (two variants): OMArk completeness vs consistency across tools.
#   A1 scatter: completeness (x) vs consistency (y), tool=colour, species=shape.
#   A2 bar    : consistency A% per tool, faceted by species (the discriminator).
# Fig B: intron precision/recall from rna_support_scoreboard.tsv.
# Parses raw OMArk .sum + scoreboard; nothing hardcoded.

suppressWarnings(suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(stringr); library(readr)
}))

repo <- normalizePath(".")
results_dir <- file.path(repo, "results")
fig_dir  <- file.path(repo, "seminar", "figures")
data_dir <- file.path(repo, "seminar", "data")
dir.create(fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

species <- c("Drosera_binata", "Drosera_paradoxa")
tools   <- c("braker_filtered", "annevo", "tiberius", "helixer")
sp_short <- function(s) str_replace(s, "Drosera_", "D. ")
tool_lab <- c(braker_filtered = "BRAKER\n(filtered)", annevo = "ANNEVO",
              tiberius = "Tiberius", helixer = "Helixer")
tool_levels <- c("braker_filtered", "annevo", "tiberius", "helixer")

sum_path <- function(sp, tool)
  file.path(results_dir, sp, "annotation", tool, "omark", "viridiplantae",
            paste0(sp, ".", tool, ".viridiplantae.sum"))

# ---- parse OMArk .sum: completeness (100-M) and consistency (A%) ----
parse_omark <- function(sp, tool) {
  f <- sum_path(sp, tool)
  if (!file.exists(f)) return(NULL)
  L <- readLines(f, warn = FALSE)
  scomp <- L[grepl("^S:[0-9.]+%", L)][1]   # S:..%,D:..%[..],M:..%
  scons <- L[grepl("^A:[0-9.]+%", L)][1]   # A:..%[..],I:..%,C:..%,U:..%
  M <- as.numeric(str_match(scomp, "M:([0-9.]+)%")[, 2])
  A <- as.numeric(str_match(scons, "^A:([0-9.]+)%")[, 2])
  U <- as.numeric(str_match(scons, "U:([0-9.]+)%")[, 2])
  tibble(species = sp, tool = tool,
         completeness = 100 - M, consistency = A, unknown = U)
}
omark <- bind_rows(lapply(species, function(sp)
            bind_rows(lapply(tools, function(t) parse_omark(sp, t))))) %>%
  mutate(sp_short = sp_short(species),
         tool = factor(tool, levels = tool_levels),
         tool_lab = tool_lab[as.character(tool)])
write_tsv(omark, file.path(data_dir, "abinitio_omark.tsv"))
message("OMArk parsed:"); print(as.data.frame(omark), digits = 4)

tool_cols <- c(braker_filtered = "#000000", annevo = "#1b9e77",
               tiberius = "#7570b3", helixer = "#d95f02")

# ---- Fig A1: scatter completeness vs consistency ----
pA1 <- ggplot(omark, aes(completeness, consistency, colour = tool, shape = sp_short)) +
  geom_point(size = 5, stroke = 1.2) +
  ggrepel::geom_text_repel(aes(label = tool_lab), size = 3.3, show.legend = FALSE,
                           segment.colour = "grey70", min.segment.length = 0) +
  scale_colour_manual(values = tool_cols, guide = "none") +
  scale_shape_manual(values = c("D. binata" = 16, "D. paradoxa" = 17), name = NULL) +
  labs(title = "OMArk: all tools are complete; they differ in consistency",
       subtitle = "Top-right is best. Completeness (conserved orthologs found) is similar across tools;\nconsistency (taxonomically coherent proteins) separates them \u2014 ANNEVO is the cleanest ab initio",
       x = "Completeness (% conserved orthologs, 100\u2212Missing)",
       y = "Consistency (% taxonomically consistent)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")
ggsave(file.path(fig_dir, "abinitio_A1_omark_scatter.png"), pA1, width = 8.5, height = 6, dpi = 200)
ggsave(file.path(fig_dir, "abinitio_A1_omark_scatter.pdf"), pA1, width = 8.5, height = 6)

# ---- Fig A2: bar, consistency + completeness per tool, faceted species ----
omk_long <- omark %>%
  select(sp_short, tool, tool_lab, Consistency = consistency, Completeness = completeness) %>%
  pivot_longer(c(Consistency, Completeness), names_to = "metric", values_to = "pct") %>%
  mutate(tool = factor(tool, levels = tool_levels),
         metric = factor(metric, levels = c("Completeness", "Consistency")))
pA2 <- ggplot(omk_long, aes(x = tool, y = pct, fill = metric)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_text(aes(label = sprintf("%.0f", pct)),
            position = position_dodge(width = 0.75), vjust = -0.3, size = 3.2) +
  facet_wrap(~ sp_short, nrow = 1) +
  scale_x_discrete(labels = tool_lab) +
  scale_fill_manual(values = c(Completeness = "#a6cee3", Consistency = "#1f78b4"), name = NULL) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25)) +
  labs(title = "OMArk completeness vs consistency by tool",
       subtitle = "Completeness is uniformly high; consistency drops for Tiberius/Helixer. ANNEVO closest to the BRAKER reference",
       x = NULL, y = "Percent") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top",
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 9))
ggsave(file.path(fig_dir, "abinitio_A2_omark_bar.png"), pA2, width = 9.5, height = 5.2, dpi = 200)
ggsave(file.path(fig_dir, "abinitio_A2_omark_bar.pdf"), pA2, width = 9.5, height = 5.2)

# ---- Fig B: intron precision/recall scatter (ab initio only) ----
ab_tools <- c("annevo", "tiberius", "helixer")
sb <- read_tsv(file.path(results_dir, "rna_support_scoreboard.tsv"), show_col_types = FALSE) %>%
  filter(tool %in% ab_tools) %>%
  mutate(sp_short = sp_short(species),
         tool = factor(tool, levels = ab_tools),
         tool_lab = tool_lab[as.character(tool)])

pB <- ggplot(sb, aes(precision, recall, colour = tool, shape = sp_short)) +
  geom_point(size = 5, stroke = 1.2) +
  ggrepel::geom_text_repel(aes(label = tool_lab), size = 3.3, show.legend = FALSE,
                           segment.colour = "grey70", min.segment.length = 0) +
  scale_colour_manual(values = tool_cols, guide = "none") +
  scale_shape_manual(values = c("D. binata" = 16, "D. paradoxa" = 17), name = NULL) +
  labs(title = "Ab initio intron support vs RNA-seq junctions",
       subtitle = "Top-right is best. Tiberius and ANNEVO are near-tied; Helixer trails on precision",
       x = "Precision (% predicted introns in RNA-seq truth set)",
       y = "Recall (% RNA-seq junctions recovered)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")
ggsave(file.path(fig_dir, "abinitio_B_intron_pr.png"), pB, width = 8.5, height = 6, dpi = 200)
ggsave(file.path(fig_dir, "abinitio_B_intron_pr.pdf"), pB, width = 8.5, height = 6)

message("done -> seminar/figures/abinitio_{A1_omark_scatter,A2_omark_bar,B_intron_pr}")
