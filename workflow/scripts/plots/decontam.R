#!/usr/bin/env Rscript
# Decontamination QC figures for the departmental seminar.
# Reads RAW per-target outputs (stats.tsv + decontamination_audit.yaml) at the
# canonical frozen-assembly targets, derives everything in-script (nothing
# hardcoded), resolves contaminant class -> coarse kingdom via NCBI taxonomy
# (taxizedb, local SQLite), and writes two figures:
#   9a  decontam_drop_dumbbell  : % contigs dropped vs % span dropped per species
#   9b  decontam_taxa_stacked   : Mb dropped per species, coloured by kingdom
#
# Canonical targets (frozen-assembly source per species):
#   aliciae, tokaiensis, binata        -> p_utg            (single)
#   scorpioides, paradoxa, roseana,
#   filiformis                         -> hap1 + hap2      (summed)
#   capensis, regia                    -> EXCLUDED (external genomes, no decontam)

suppressPackageStartupMessages({
  library(yaml); library(dplyr); library(tidyr)
  library(ggplot2); library(stringr); library(purrr); library(readr)
})

repo <- normalizePath(".")
results_dir <- file.path(repo, "results")
fig_dir  <- file.path(repo, "seminar", "figures")
data_dir <- file.path(repo, "seminar", "data")
dir.create(fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# ---- canonical targets ----------------------------------------------------
targets <- tribble(
  ~species,               ~haps,
  "Drosera_aliciae",      "p_utg",
  "Drosera_tokaiensis",   "p_utg",
  "Drosera_binata",       "p_utg",
  "Drosera_scorpioides",  "hap1,hap2",
  "Drosera_paradoxa",     "hap1,hap2",
  "Drosera_roseana",      "hap1,hap2",
  "Drosera_filiformis",   "hap1,hap2",
) %>% separate_rows(haps, sep = ",")

stats_path <- function(sp, h)
  file.path(results_dir, sp, "assembly", "decontaminated", h, "stats.tsv")
audit_path <- function(sp, h)
  file.path(results_dir, sp, "assembly", "decontaminated", h, "decontamination_audit.yaml")

# ---- 9a: read + sum stats.tsv per species ---------------------------------
read_stats <- function(p) {
  x <- readr::read_tsv(p, show_col_types = FALSE)
  setNames(as.numeric(x$value), x$metric)
}

stats9a <- targets %>%
  mutate(s = map2(species, haps, ~ read_stats(stats_path(.x, .y)))) %>%
  mutate(total_contigs = map_dbl(s, "total_contigs"),
         dropped_contigs = map_dbl(s, "dropped_contigs"),
         total_span = map_dbl(s, "total_span"),
         dropped_span = map_dbl(s, "dropped_span")) %>%
  group_by(species) %>%
  summarise(total_contigs = sum(total_contigs),
            dropped_contigs = sum(dropped_contigs),
            total_span = sum(total_span),
            dropped_span = sum(dropped_span), .groups = "drop") %>%
  mutate(drop_contig_pct = 100 * dropped_contigs / total_contigs,
         drop_span_pct   = 100 * dropped_span   / total_span,
         dropped_span_mb = dropped_span / 1e6,
         sp_short = str_replace(species, "Drosera_", "D. "))

write_tsv(stats9a, file.path(data_dir, "decontam_9a.tsv"))
message("9a table:"); print(as.data.frame(stats9a[, c("species","drop_contig_pct","drop_span_pct","dropped_span_mb")]), digits = 3)

# ---- 9b: parse audit.yaml drop_class_span per species ---------------------
read_drop_classes <- function(sp, h) {
  y <- yaml::read_yaml(audit_path(sp, h))
  dcs <- y$results$drop_class_span
  if (is.null(dcs) || length(dcs) == 0) return(tibble(class = character(), bp = numeric()))
  tibble(class = names(dcs), bp = as.numeric(unlist(dcs)))
}

drops_long <- targets %>%
  mutate(d = map2(species, haps, read_drop_classes)) %>%
  select(species, d) %>% unnest(d) %>%
  group_by(species, class) %>% summarise(bp = sum(bp), .groups = "drop")

# ---- taxonomy: class -> coarse kingdom via taxizedb (local NCBI db) -------
suppressPackageStartupMessages(library(taxizedb))
# normalise blobtools class names: strip "-undef" suffixes, keep base name
clean_class <- function(x) str_replace(x, "-undef$", "")
uniq_classes <- unique(clean_class(drops_long$class))

message(sprintf("resolving %d unique classes via taxizedb...", length(uniq_classes)))
taxids <- taxizedb::name2taxid(uniq_classes, out_type = "summary")  # may return multiple
# name2taxid summary: data.frame(name, id). Keep first id per name.
taxids <- taxids %>% group_by(name) %>% slice(1) %>% ungroup()

cls <- taxizedb::classification(taxids$id, db = "ncbi")
# coarse kingdom from lineage: prefer kingdom rank; bacteria/archaea -> superkingdom
coarse_kingdom <- function(lin) {
  if (is.null(lin) || !is.data.frame(lin) || all(is.na(lin$name))) return(NA_character_)
  nm <- lin$name
  if (any(nm %in% c("Bacteria","Pseudomonadati","Bacillati","Thermotogati",
                    "Fusobacteriati","Thermodesulfobacteriati","Nitrospirati",
                    "Spirochaetati","Fibrobacterati","Planctomycetati"))) return("Bacteria")
  if (any(nm == "Archaea")) return("Archaea")
  if (any(nm %in% c("Viruses","Orthornavirae","Bamfordvirae","Heunggongvirae",
                    "Shotokuvirae","Pararnavirae"))) return("Virus")
  kg <- lin$name[lin$rank == "kingdom"]
  if (length(kg)) {
    if (kg %in% "Fungi")         return("Fungi")
    if (kg %in% "Metazoa")       return("Metazoa")
    if (kg %in% "Viridiplantae") return("Viridiplantae")
  }
  "Other-Eukaryote"
}
king_map <- tibble(
  id = names(cls),
  kingdom = map_chr(cls, coarse_kingdom)
) %>% left_join(taxids, by = "id") %>%
  transmute(class = name, kingdom)

# write the lookup for human inspection (you asked to be able to eyeball it)
write_tsv(king_map, file.path(data_dir, "class_kingdom_map.tsv"))
message("class->kingdom map (inspect seminar/data/class_kingdom_map.tsv):")
print(as.data.frame(king_map))

# ---- join + per-species kingdom Mb ----------------------------------------
king_levels <- c("Fungi","Metazoa","Viridiplantae","Bacteria","Other-Eukaryote","Archaea","Virus","Other")
king_cols   <- c(Fungi="#b15928", Metazoa="#1f78b4", Viridiplantae="#33a02c",
                 Bacteria="#6a3d9a", `Other-Eukaryote`="#ff7f00",
                 Archaea="#a6cee3", Virus="#e31a1c", Other="#999999")

drops_king <- drops_long %>%
  mutate(class = clean_class(class)) %>%
  left_join(king_map, by = "class") %>%
  mutate(kingdom = ifelse(is.na(kingdom), "Other", kingdom),
         kingdom = factor(kingdom, levels = king_levels),
         mb = bp / 1e6,
         sp_short = str_replace(species, "Drosera_", "D. ")) %>%
  group_by(species, sp_short, kingdom) %>%
  summarise(mb = sum(mb), .groups = "drop")

write_tsv(drops_king, file.path(data_dir, "decontam_9b_kingdom.tsv"))

# order species by total dropped Mb (descending) for both plots
sp_order <- stats9a %>% arrange(desc(dropped_span_mb)) %>% pull(sp_short)

# ---- PLOT 9a: dumbbell ----------------------------------------------------
d9a <- stats9a %>%
  mutate(sp_short = factor(sp_short, levels = rev(sp_order))) %>%
  select(sp_short, `% contigs` = drop_contig_pct, `% span` = drop_span_pct)

d9a_long <- d9a %>% pivot_longer(-sp_short, names_to = "metric", values_to = "pct")

p9a <- ggplot(d9a, aes(y = sp_short)) +
  geom_segment(aes(x = `% span`, xend = `% contigs`, yend = sp_short),
               colour = "grey70", linewidth = 1.2) +
  geom_point(data = d9a_long, aes(x = pct, colour = metric), size = 4) +
  scale_colour_manual(values = c("% contigs" = "#d95f02", "% span" = "#1b9e77"),
                      name = NULL) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(title = "Decontamination: many contigs removed, little sequence lost",
       subtitle = "Dropped contigs are short fragments \u2014 the bulk of the assembly is clean",
       x = "% of assembly dropped", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "decontam_9a_dropped.png"), p9a, width = 8, height = 4.5, dpi = 200)
ggsave(file.path(fig_dir, "decontam_9a_dropped.pdf"), p9a, width = 8, height = 4.5)

# ---- PLOT 9b: stacked kingdom Mb ------------------------------------------
d9b <- drops_king %>% mutate(sp_short = factor(sp_short, levels = rev(sp_order)))

p9b <- ggplot(d9b, aes(x = mb, y = sp_short, fill = kingdom)) +
  geom_col() +
  scale_fill_manual(values = king_cols, name = "Contaminant\nkingdom", drop = TRUE) +
  labs(title = "Contaminant composition across the genus",
       subtitle = "Recurring Fungi & Arthropoda (Metazoa); species-specific dominance varies",
       x = "Dropped sequence (Mb)", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))



# ---- PLOT 9c: recurring taxa across species ----
# ---- PLOT 9c: recurring taxa across species ----
recurring <- drops_long %>%
  mutate(class = clean_class(class)) %>%
  left_join(king_map, by = "class") %>%
  mutate(kingdom = ifelse(is.na(kingdom), "Other", kingdom),
         mb = bp / 1e6) %>%
  group_by(class, kingdom) %>%
  summarise(species_count = n_distinct(species), total_mb = sum(mb), .groups = "drop") %>%
  arrange(desc(total_mb)) %>% slice_head(n = 15) %>%
  mutate(class = factor(class, levels = rev(class)))


readr::write_tsv(recurring, file.path(data_dir, "decontam_9c_recurring.tsv"))
p9c <- ggplot(recurring, aes(x = total_mb, y = class, fill = kingdom)) +
  geom_col() +
  geom_text(aes(label = sprintf("%d sp", species_count)), hjust = -0.15, size = 3.4, colour = "grey30") +
  scale_fill_manual(values = king_cols, name = "Kingdom", drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Recurring contaminant taxa across the genus",
       subtitle = "Top dropped classes by total span; label = number of species affected",
       x = "Total dropped sequence across species (Mb)", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold"))
ggsave(file.path(fig_dir, "decontam_9c_recurring.png"), p9c, width = 8.5, height = 5, dpi = 200)
ggsave(file.path(fig_dir, "decontam_9c_recurring.pdf"), p9c, width = 8.5, height = 5)


ggsave(file.path(fig_dir, "decontam_9b_taxa.png"), p9b, width = 8.5, height = 4.5, dpi = 200)
ggsave(file.path(fig_dir, "decontam_9b_taxa.pdf"), p9b, width = 8.5, height = 4.5)

message("done. figures -> seminar/figures/  data -> seminar/data/")
