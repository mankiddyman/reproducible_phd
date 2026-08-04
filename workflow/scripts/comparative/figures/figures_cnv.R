#!/usr/bin/env Rscript
# =============================================================================
# figures_cnv.R  --  presentation figures for the holocentricity CNV work
#   Fig 1: RAD21/SCC1 cohesin copy number, holo vs mono (+ locus annotation)
#   Fig 2: ploidy-inference validation (mode recapitulates known assembly ploidy)
# Run from repo root. Writes PDF + PNG to results/comparative/holocentricity/figures/
# =============================================================================

suppressWarnings(suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
}))

repo    <- normalizePath(".")
gc_file <- file.path(repo, "results/comparative/orthofinder/out/Results_drosera/Orthogroups/Orthogroups.GeneCount.tsv")
fig_dir <- file.path(repo, "results/comparative/holocentricity/figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# species order (holo first, then mono) + short labels + type
sp_levels <- c("Drosera_paradoxa","Drosera_regia","Drosera_roseana",
               "Drosera_aliciae","Drosera_binata","Drosera_capensis",
               "Drosera_filiformis","Drosera_tokaiensis")
sp_short  <- c("paradoxa","regia","roseana","aliciae","binata","capensis","filiformis","tokaiensis")
sp_type   <- c("holocentric","holocentric","holocentric",
               "monocentric","monocentric","monocentric","monocentric","monocentric")
names(sp_short) <- sp_levels; names(sp_type) <- sp_levels

# shared theme
th <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top",
        plot.title = element_text(face = "plain", size = 14),
        plot.subtitle = element_text(size = 11, colour = "grey40"))
col_type <- c(holocentric = "#1D9E75", monocentric = "#BA7517")  # teal / amber

# =============================================================================
# FIG 1 — RAD21/SCC1 cohesin copy number per species
# =============================================================================
# RAW member counts for OG0001569 (from Orthogroups.tsv), with the locus
# interpretation we verified by scaffold (holo "2" = hap1/hap2 of ONE locus).
rad21 <- data.table(
  species   = sp_levels,
  raw_copies = c(2, 1, 2,    6, 4, 4, 3, 6),         # raw OrthoFinder member counts
  loci       = c(1, 1, 1,    NA, 2, 2, NA, NA)       # verified distinct loci (NA = not scaffold-checked / polyploid)
)
rad21[, short := factor(sp_short[species], levels = sp_short)]
rad21[, type  := sp_type[species]]

p1 <- ggplot(rad21, aes(short, raw_copies, fill = type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = raw_copies), vjust = -0.4, size = 4) +
  scale_fill_manual(values = col_type, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "RAD21/SCC1 cohesin kleisin — copy number across Drosera",
       subtitle = "holocentrics retain 1 locus (phased = hap1/hap2 of one gene); monocentrics carry >=2 loci on different chromosomes",
       x = NULL, y = "gene copies (raw)") +
  th
ggsave(file.path(fig_dir, "fig1_rad21_copies.pdf"), p1, width = 7.5, height = 4.5)
ggsave(file.path(fig_dir, "fig1_rad21_copies.png"), p1, width = 7.5, height = 4.5, dpi = 200)
cat("wrote fig1_rad21_copies.{pdf,png}\n")

# =============================================================================
# FIG 2 — ploidy-inference validation
# =============================================================================
# load gene counts, restrict to OGs present in all 8, compute the per-species
# copy-number distribution + inferred mode (this is exactly what cnv_ploidy.R does).
gc <- fread(gc_file); if ("Total" %in% names(gc)) gc[, Total := NULL]
gc <- gc[, c("Orthogroup", sp_levels), with = FALSE]
present_all <- gc[, rowSums(.SD > 0) == length(sp_levels), .SDcols = sp_levels]
cnv <- gc[present_all == TRUE]

# long format of copy numbers (cap at 8 for display)
long <- melt(cnv, id.vars = "Orthogroup", measure.vars = sp_levels,
             variable.name = "species", value.name = "copies")
long[, copies_disp := pmin(copies, 8)]
long[, short := factor(sp_short[as.character(species)], levels = sp_short)]
long[, type  := sp_type[as.character(species)]]

# inferred mode + known ploidy for annotation
mode_dt <- long[copies >= 1, .(mode = as.integer(names(which.max(table(copies))))), by = short]
known   <- data.table(short = factor(sp_short, levels = sp_short),
                       known = c(2,1,2, 3,2,2,2,3))   # paradoxa2 regia1 roseana2 | ali3 bin2 cap2 fil2 tok3
ann <- merge(mode_dt, known, by = "short")
ann[, lab := paste0("inferred ", mode, " / known ", known,
                    ifelse(mode == known, "  match", "  MISMATCH"))]

# small-multiples histogram of copy number per species, mode line
p2 <- ggplot(long[copies >= 1], aes(copies_disp, fill = type)) +
  geom_bar(width = 0.9) +
  geom_vline(data = mode_dt, aes(xintercept = mode), colour = "red", linewidth = 0.5, linetype = "22") +
  geom_text(data = ann, aes(x = 5.2, y = Inf, label = lab),
            inherit.aes = FALSE, vjust = 1.6, size = 3, colour = "grey25") +
  facet_wrap(~ short, ncol = 4) +
  scale_fill_manual(values = col_type, name = NULL) +
  scale_x_continuous(breaks = 1:8, labels = c(1:7, "8+")) +
  labs(title = "Empirical ploidy inference recapitulates known assembly ploidy",
       subtitle = "copy-number distribution per species (OGs present in all 8); red line = modal copy number = inferred baseline ploidy",
       x = "copies per orthogroup", y = "orthogroups") +
  th + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
ggsave(file.path(fig_dir, "fig2_ploidy_validation.pdf"), p2, width = 9, height = 5)
ggsave(file.path(fig_dir, "fig2_ploidy_validation.png"), p2, width = 9, height = 5, dpi = 200)
cat("wrote fig2_ploidy_validation.{pdf,png}\n")

cat("\nfigures in:", fig_dir, "\n")
