#!/usr/bin/env Rscript
# Slide 10.6 — OMArk across species: decontamination leaves 0% contamination
# everywhere; the remaining differentiator is assembly type (dualhap non-HiC
# assemblies carry more "unknown" / unplaceable proteins).
# A + I + C + U sum to 100 (whole-proteome denominator) -> a stacked bar is valid.
# Variant (a): all species, annotation method labelled (braker vs annevo).
# Variant (b): annevo-annotated non-HiC four only (collapsed vs dualhap) -
#              method controlled, so the assembly-type effect is unconfounded.

suppressWarnings(suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(stringr); library(readr); library(purrr)
}))

repo <- normalizePath(".")
results_dir <- file.path(repo, "results")
fig_dir  <- file.path(repo, "seminar", "figures")
data_dir <- file.path(repo, "seminar", "data")
dir.create(fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# species -> (annotation method actually used, assembly class)
meta <- tribble(
  ~species,              ~tool,             ~assembly,
  "Drosera_binata",      "braker_filtered", "Chromosome level + (RNA)",
  "Drosera_paradoxa",    "braker_filtered", "Chromosome level + (RNA)",
  "Drosera_capensis",    "braker_filtered", "Chromosome level + (RNA)",
  "Drosera_regia",       "braker_filtered", "Chromosome level + (RNA)",
  "Drosera_aliciae",     "annevo",          "Contig level polyploid",
  "Drosera_tokaiensis",  "annevo",          "Contig level polyploid",
  "Drosera_roseana",     "annevo",          "Contig level diploid",
  "Drosera_filiformis",  "annevo",          "Contig level diploid"
)
sp_short <- function(s) str_replace(s, "Drosera_", "D. ")
sum_path <- function(sp, tool)
  file.path(results_dir, sp, "annotation", tool, "omark", "viridiplantae",
            paste0(sp, ".", tool, ".viridiplantae.sum"))

parse_cons <- function(sp, tool) {
  f <- sum_path(sp, tool); if (!file.exists(f)) return(NULL)
  L <- readLines(f, warn = FALSE)
  l  <- L[grepl("^A:[0-9.]+%", L)][1]   # A,I,C,U  (whole proteome)
  lc <- L[grepl("^S:[0-9.]+%", L)][1]   # S,D,M    (conserved HOGs)
  num <- function(line, tag) as.numeric(str_match(line, paste0(tag, ":([0-9.]+)%"))[, 2])
  tibble(Consistent = num(l,"^A"), Inconsistent = num(l,"I"),
         Contamination = num(l,"C"), Unknown = num(l,"U"),
         completeness = 100 - num(lc,"M"), consistency = num(l,"^A"))
}

omk <- meta %>%
  mutate(d = map2(species, tool, parse_cons)) %>%
  select(species, tool, assembly, d) %>% unnest(d) %>%
  mutate(sp_short = sp_short(species))
write_tsv(omk, file.path(data_dir, "omark_noHiC.tsv"))
message("OMArk consistency table:"); print(as.data.frame(omk), digits = 4)
message(sprintf("max contamination across all species: %.2f%%", max(omk$Contamination)))

cat_levels <- c("Consistent", "Inconsistent", "Unknown")  # C omitted (0 everywhere)
cat_cols   <- c(Consistent = "#1b9e77", Inconsistent = "#999999", Unknown = "#d95f02")

to_long <- function(df) df %>%
  select(sp_short, assembly, tool, all_of(cat_levels)) %>%
  pivot_longer(all_of(cat_levels), names_to = "cat", values_to = "pct") %>%
  mutate(cat = factor(cat, levels = rev(cat_levels)))

# assembly ordering for the y-axis (group related types together)
asm_order_a <- c("Chromosome level + (RNA)", "Contig level polyploid", "Contig level diploid")

# ---- Variant (a): all species, method labelled ----
da <- omk %>% mutate(assembly = factor(assembly, levels = asm_order_a)) %>%
  arrange(assembly, Unknown) %>%
  mutate(lab = paste0(sp_short, "  [", ifelse(tool == "braker_filtered", "BRAKER", "ANNEVO"), "]"),
         lab = factor(lab, levels = rev(unique(lab))))
da_long <- da %>% select(lab, assembly, all_of(cat_levels)) %>%
  pivot_longer(all_of(cat_levels), names_to = "cat", values_to = "pct") %>%
  mutate(cat = factor(cat, levels = rev(cat_levels)))

pa <- ggplot(da_long, aes(x = pct, y = lab, fill = cat)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = cat_cols, name = NULL, breaks = cat_levels) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(title = "OMArk proteome consistency across all species",
       subtitle = "Contamination = 0% everywhere; unknown proteins rise in non-HiC dualhap (BRAKER vs ANNEVO confounded)",
       x = "% of predicted proteome", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top",
        panel.grid.major.y = element_blank())
ggsave(file.path(fig_dir, "omark_noHiC_a_all.png"), pa, width = 9.5, height = 5, dpi = 200)
ggsave(file.path(fig_dir, "omark_noHiC_a_all.pdf"), pa, width = 9.5, height = 5)

# ---- Variant (b): annevo-only non-HiC four (method controlled) ----
db <- omk %>% filter(tool == "annevo") %>%
  mutate(assembly = factor(assembly, levels = c("Contig level polyploid", "Contig level diploid"))) %>%
  arrange(assembly, Unknown) %>%
  mutate(lab = paste0(sp_short, "  [", str_replace(assembly, "Contig level ", ""), "]"),
         lab = factor(lab, levels = rev(unique(lab))))
db_long <- db %>% select(lab, all_of(cat_levels)) %>%
  pivot_longer(all_of(cat_levels), names_to = "cat", values_to = "pct") %>%
  mutate(cat = factor(cat, levels = rev(cat_levels)))

pb <- ggplot(db_long, aes(x = pct, y = lab, fill = cat)) +
  geom_col(width = 0.66) +
  scale_fill_manual(values = cat_cols, name = NULL, breaks = cat_levels) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(title = "Non-HiC annotations (ANNEVO): dualhap assemblies are noisier",
       subtitle = "Same method, all decontaminated. Dualhap (roseana, filiformis) ~25% unplaceable vs ~11% collapsed",
       x = "% of predicted proteome", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top",
        panel.grid.major.y = element_blank())
ggsave(file.path(fig_dir, "omark_noHiC_b_annevo.png"), pb, width = 9, height = 3.8, dpi = 200)
ggsave(file.path(fig_dir, "omark_noHiC_b_annevo.pdf"), pb, width = 9, height = 3.8)

# ---- Variant (c): completeness vs consistency scatter, all species ----
# For binata & paradoxa we ALSO plot their ANNEVO point and draw an arrow
# ANNEVO -> BRAKER, visualising what RNA-informed annotation buys you.
library(ggrepel)
asm_cols <- c("Chromosome level + (RNA)" = "#1b9e77",
              "Contig level polyploid" = "#e6ab02",
              "Contig level diploid" = "#d95f02")

# base points: each species' primary (frozen) annotation, as before
sc <- omk %>%
  mutate(assembly = factor(assembly, levels = names(asm_cols)),
         method = ifelse(tool == "braker_filtered", "BRAKER", "ANNEVO"),
         lab = sp_short(species))

# extra: ANNEVO points for the two RNA species (so the arrow has a tail)
extra_sp <- c("Drosera_binata", "Drosera_paradoxa")
annevo_pts <- bind_rows(lapply(extra_sp, function(spp) {
  d <- parse_cons(spp, "annevo")
  d$species <- spp; d
})) %>%
  mutate(assembly = "Chromosome level + (RNA)",  # same genome, just ab-initio annotation
         method = "ANNEVO", lab = sp_short(species),
         assembly = factor(assembly, levels = names(asm_cols)))

# arrow segments: ANNEVO (tail) -> BRAKER (head) per RNA species
braker_pts <- sc %>% filter(species %in% extra_sp) %>%
  select(species, bx = completeness, by = consistency)
arrows <- annevo_pts %>% select(species, ax = completeness, ay = consistency) %>%
  left_join(braker_pts, by = "species")

all_pts <- bind_rows(sc, annevo_pts)

pc <- ggplot(all_pts, aes(completeness, consistency, colour = assembly, shape = method)) +
  geom_segment(data = arrows, inherit.aes = FALSE,
               aes(x = ax, y = ay, xend = bx, yend = by),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               colour = "grey50", linewidth = 0.6) +
  geom_point(size = 5, stroke = 1.3) +
  geom_text_repel(data = sc, aes(label = lab), size = 3.4, show.legend = FALSE,
                  segment.colour = "grey70", min.segment.length = 0, max.overlaps = Inf) +
  scale_colour_manual(values = asm_cols, name = "Assembly") +
  scale_shape_manual(values = c("BRAKER" = 16, "ANNEVO" = 17), name = "Method") +
  labs(title = "Genome quality space: assembly tier and RNA both shape quality",
       subtitle = "Top-right is best. Arrows: ANNEVO -> BRAKER for the two RNA species (RNA lifts consistency). Contamination = 0% in all",
       x = "Completeness (% conserved orthologs)",
       y = "Consistency (% taxonomically consistent)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")
ggsave(file.path(fig_dir, "omark_noHiC_c_scatter.png"), pc, width = 9, height = 5.5, dpi = 200)
ggsave(file.path(fig_dir, "omark_noHiC_c_scatter.pdf"), pc, width = 9, height = 5.5)

message("done -> seminar/figures/omark_noHiC_{a_all,b_annevo}.{png,pdf}")
