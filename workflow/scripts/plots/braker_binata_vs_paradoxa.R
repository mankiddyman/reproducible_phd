#!/usr/bin/env Rscript
# Slide 10.1 — Why big genomes need braker_filtering (binata vs paradoxa).
# Figure 1: raw-BRAKER support metrics (single-exon %, intron support %, exon
#           support %) per species — the problem.
# Figure 2: raw -> filtered, transcripts split single/multi-exon — the fix
#           (filtering removes the single-exon junk; binata barely changes).
# Everything derived from raw files: gene_support.tsv (raw braker per-transcript
# table) + the filtered omark_clean gff3 (survivor transcript IDs). No hardcoding.

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
sp_short <- function(s) str_replace(s, "Drosera_", "D. ")

gs_path <- function(sp) file.path(results_dir, sp, "annotation", "braker",
                                  "output", sp, "results", "gene_support.tsv")
filt_gff <- function(sp) file.path(results_dir, sp, "annotation", "braker_filtered",
                                   paste0(sp, ".braker_filtered.omark_clean.gff3"))

# ---- read raw braker per-transcript support table (skip # comment lines) ----
read_gs <- function(sp) {
  readr::read_tsv(gs_path(sp), comment = "#", show_col_types = FALSE,
                  na = c("NA", "")) %>%
    mutate(species = sp)
}
gs <- bind_rows(lapply(species, read_gs))

# ---- FIGURE 1: raw-braker support metrics, 3 metrics x 2 species ----
metrics <- gs %>% group_by(species) %>%
  summarise(
    n_tx          = n(),
    single_exon_pct = 100 * mean(num_introns == 0),
    intron_sup_pct  = 100 * sum(introns_sup_any, na.rm = TRUE) / sum(num_introns, na.rm = TRUE),
    exon_sup_pct    = 100 * sum(exons_sup_any,  na.rm = TRUE) / sum(num_exons,  na.rm = TRUE),
    .groups = "drop") %>%
  mutate(sp_short = sp_short(species))

write_tsv(metrics, file.path(data_dir, "braker_101_metrics.tsv"))
message("Figure 1 metrics:"); print(as.data.frame(metrics), digits = 4)

m_long <- metrics %>%
  select(sp_short,
         `Single-exon\ntranscripts` = single_exon_pct,
         `Intron support\n(of all introns)` = intron_sup_pct,
         `Exon support\n(of all exons)` = exon_sup_pct) %>%
  pivot_longer(-sp_short, names_to = "metric", values_to = "pct") %>%
  mutate(metric = factor(metric, levels = c("Single-exon\ntranscripts",
                                            "Intron support\n(of all introns)",
                                            "Exon support\n(of all exons)")))

p1 <- ggplot(m_long, aes(x = metric, y = pct, fill = sp_short)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            position = position_dodge(width = 0.75), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("D. binata" = "#1b9e77", "D. paradoxa" = "#d95f02"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25)) +
  labs(title = "Raw BRAKER predictions: the big genome is noisier",
       subtitle = "D. paradoxa (4.7 Gb) predicts far more single-exon genes; introns are equally well-supported in both,\nbut paradoxa's extra single-exon models drag down exon support",
       x = NULL, y = "Percent") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top",
        panel.grid.major.x = element_blank())

ggsave(file.path(fig_dir, "braker_101_metrics.png"), p1, width = 9, height = 5.2, dpi = 200)
ggsave(file.path(fig_dir, "braker_101_metrics.pdf"), p1, width = 9, height = 5.2)

# ---- FIGURE 2: raw -> filtered, single/multi-exon composition ----
survivors <- function(sp) {
  gff <- filt_gff(sp)
  ln <- readLines(gff, warn = FALSE)
  ln <- ln[!startsWith(ln, "#")]
  tx <- ln[grepl("\t(transcript|mRNA)\t", ln)]
  ids <- str_match(tx, "ID=([^;]+)")[, 2]
  unique(ids)
}

comp <- bind_rows(lapply(species, function(sp) {
  surv <- survivors(sp)
  g <- gs %>% filter(species == sp) %>%
    mutate(exon_class = ifelse(num_introns == 0, "Single-exon", "Multi-exon"),
           kept = transcript_id %in% surv)
  bind_rows(
    g %>% mutate(stage = "Raw BRAKER"),
    g %>% filter(kept) %>% mutate(stage = "braker_filtered")
  ) %>%
    count(species, stage, exon_class) %>%
    mutate(sp_short = sp_short(species))
}))

dropped <- bind_rows(lapply(species, function(sp) {
  surv <- survivors(sp)
  gs %>% filter(species == sp) %>%
    mutate(exon_class = ifelse(num_introns == 0, "Single-exon", "Multi-exon"),
           dropped = !(transcript_id %in% surv)) %>%
    filter(dropped) %>% count(species, exon_class)
}))
dropped_summary <- dropped %>% group_by(species) %>%
  mutate(total = sum(n), pct = 100 * n / total) %>% ungroup()
message("dropped-transcript composition:"); print(as.data.frame(dropped_summary), digits = 3)
write_tsv(comp, file.path(data_dir, "braker_101_composition.tsv"))

comp <- comp %>%
  mutate(stage = factor(stage, levels = c("Raw BRAKER", "braker_filtered")),
         exon_class = factor(exon_class, levels = c("Multi-exon", "Single-exon")))

p2 <- ggplot(comp, aes(x = stage, y = n, fill = exon_class)) +
  geom_col(width = 0.7) +
  facet_wrap(~ sp_short, nrow = 1) +
  scale_fill_manual(values = c("Single-exon" = "#d95f02", "Multi-exon" = "#7570b3"),
                    name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Filtering to RNA-supported models removes single-exon junk",
       subtitle = "D. paradoxa sheds tens of thousands of single-exon transcripts; D. binata was already clean",
       x = NULL, y = "Transcripts") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top",
        panel.grid.major.x = element_blank())

ggsave(file.path(fig_dir, "braker_101_filter.png"), p2, width = 8.5, height = 5, dpi = 200)
ggsave(file.path(fig_dir, "braker_101_filter.pdf"), p2, width = 8.5, height = 5)

message("done. figures -> seminar/figures/braker_101_{metrics,filter}.{png,pdf}")
