#!/usr/bin/env Rscript
# =====================================================================
# Presence/Absence holocentricity screen
# ---------------------------------------------------------------------
# Question: which orthogroups are GAINED or LOST in association with
# holocentricity (present in one centromere group, absent in the other)?
#
# Why P/A and not copy number: presence (>=1 gene) is robust to the
# annotation-DEPTH / assembly-tier inflation that confounds raw copy
# number (e.g. roseana's gene-rich dualhap assembly). It is NOT fully
# robust to annotation INCOMPLETENESS: a lean annotation (regia ~22k,
# capensis ~25k genes) can show false absences. OMArk completeness ~90%
# bounds this, but treat absences in lean annotations with caution.
#
# Stats: per orthogroup, 2x2 (present/absent x holo/mono), Fisher exact.
# n = 3 holo (paradoxa, roseana, regia) vs 5 mono. With only 8 species
# there are 4x6=24 possible presence patterns -> 24 distinct p-values;
# min achievable p ~0.018 (one extreme pattern). These are CANDIDATES,
# not biomarkers. scorpioides excluded (genome pending).
# =====================================================================
suppressWarnings(suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
}))

repo <- normalizePath(".")
gc_file <- file.path(repo, "results/comparative/orthofinder/out/Results_drosera/Orthogroups/Orthogroups.GeneCount.tsv")
out_dir <- file.path(repo, "results/comparative/holocentricity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- labels ----
# read centromere from species.csv where available; capensis/regia are
# external genomes (not in species.csv) -> set explicitly.
sp_csv <- fread(file.path(repo, "config/species.csv"))
cent <- setNames(sp_csv$centromere, sp_csv$species_id)
cent["Drosera_capensis"] <- "monocentric"   # external genome
cent["Drosera_regia"]    <- "holocentric"   # external genome

# ---- copy-number matrix ----
gc <- fread(gc_file)
gc[, Total := NULL]
sp_cols <- setdiff(names(gc), "Orthogroup")

# align labels to matrix columns
stopifnot(all(sp_cols %in% names(cent)))
holo_sp <- sp_cols[cent[sp_cols] == "holocentric"]
mono_sp <- sp_cols[cent[sp_cols] == "monocentric"]
n_holo <- length(holo_sp); n_mono <- length(mono_sp)
cat("Holo (", n_holo, "):", paste(holo_sp, collapse=", "), "\n")
cat("Mono (", n_mono, "):", paste(mono_sp, collapse=", "), "\n\n")

# ---- presence counts per orthogroup ----
pres <- gc[, lapply(.SD, function(x) as.integer(x > 0)), .SDcols = sp_cols]
dt <- data.table(
  Orthogroup = gc$Orthogroup,
  a = rowSums(pres[, ..holo_sp]),   # holo species present
  b = rowSums(pres[, ..mono_sp])    # mono species present
)
dt[, n_present := a + b]

# ---- Fisher per unique (a,b) pattern ----
pat <- unique(dt[, .(a, b)])
pat[, pid := .I]
fish <- pat[, {
  m <- matrix(c(a, n_holo - a, b, n_mono - b), nrow = 2)  # rows present/absent, cols holo/mono
  ft <- suppressWarnings(fisher.test(m))
  .(p_fisher = ft$p.value, odds_ratio = unname(ft$estimate))
}, by = pid]
pat <- merge(pat, fish, by = "pid")
pat[, delta := a / n_holo - b / n_mono]        # presence-fraction diff, -1..1
pat[, direction := fifelse(delta > 0, "holo-enriched",
                    fifelse(delta < 0, "mono-enriched", "no difference"))]

dt <- merge(dt, pat[, .(a, b, p_fisher, odds_ratio, delta, direction)],
            by = c("a", "b"))
dt[, q_bh := p.adjust(p_fisher, method = "BH")]
dt[, holo_frac := a / n_holo][, mono_frac := b / n_mono]
dt[, abs_delta := abs(delta)]
setorder(dt, p_fisher, -abs_delta)

fwrite(dt, file.path(out_dir, "pa_screen.tsv"), sep = "\t")

# ---- summary ----
cat("Total orthogroups:", nrow(dt), "\n")
cat("Present in all 8 (uninformative core):", dt[a==n_holo & b==n_mono, .N], "\n")
cat("p_fisher <= 0.05 :", dt[p_fisher <= 0.05, .N], "\n")
cat("q_bh    <= 0.10 :", dt[q_bh <= 0.10, .N], "\n\n")

cat("=== presence-pattern table (the 24 cells; n = orthogroups) ===\n")
patt_summary <- dt[, .N, by = .(a, b, holo_frac, mono_frac, delta, p_fisher, q_bh, direction)][order(p_fisher, -abs(delta))]  # base order() allows abs()
print(patt_summary[, .(holo=paste0(a,"/",n_holo), mono=paste0(b,"/",n_mono),
                       delta=round(delta,2), p=round(p_fisher,4),
                       q=round(q_bh,3), direction, n_OGs=N)])


cat("\n=== all or nothing candidates (present in all holo, absent in all mono, or vice versa) ===\n")
print(dt[(delta == 1 | delta == -1), .(Orthogroup, a, b, delta=round(delta,2), p_fisher=round(p_fisher,4), q_bh=round(q_bh,3))])
print(paste(" all holo candidates:", dt[delta == 1, .N], "(", round(dt[delta == 1, .N] / nrow(dt) * 100, 2), "% of all orthogroups)"))
print(paste(" all mono candidates:", dt[delta == -1, .N], "(", round(dt[delta == -1, .N] / nrow(dt) * 100, 2), "% of all orthogroups)"))

# write only those collumns to a file called pa_complete_candidates.tsv
fwrite(dt[(a==n_holo & b==0) | (a==0 & b==n_mono), .(Orthogroup, a, b, delta=round(delta,2), p_fisher=round(p_fisher,4), q_bh=round(q_bh,3))], file.path(out_dir, "pa_complete_candidates.tsv"), sep = "\t")


# ---- volcano: one point per pattern, sized by n orthogroups ----
vol <- dt[, .N, by = .(delta, p_fisher, direction, a, b)]
vol[, neglog10p := -log10(p_fisher)]
vol[, lab := paste0(a, "h/", b, "m")]
vol[, candidate := (a == n_holo & b == 0) | (a == 0 & b == n_mono)]

n_total <- nrow(dt)
vol[, pct := round(N / n_total * 100, 2)]
# label text: pattern + count + %  (e.g. "3h/0m\n105 OGs (0.38%)")
vol[, lab_full := sprintf("%s\n%d OGs (%.2f%%)", lab, N, pct)]

dir_cols <- c("holo-enriched" = "#1b9e77", "mono-enriched" = "#d95f02", "no difference" = "#999999")
vol[, fill_col := ifelse(candidate, dir_cols[direction], NA)]

p <- ggplot(vol, aes(delta, neglog10p, size = N, colour = direction)) +
  geom_point(shape = 21, stroke = 1.1, aes(fill = I(fill_col))) +
  ggrepel::geom_text_repel(data = vol[candidate == TRUE], aes(label = lab_full),
                           size = 3.3, fontface = "bold", show.legend = FALSE,
                           lineheight = 0.9, segment.colour = "grey70",
                           max.overlaps = Inf, point.padding = 0.5, box.padding = 0.6) +
  scale_colour_manual(values = dir_cols, name = NULL) +
  scale_size_continuous(range = c(2, 12), name = "n orthogroups") +
  scale_x_continuous(limits = c(-1, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  labs(title = "Presence/absence association with holocentricity",
       subtitle = sprintf("%d holo vs %d mono species; one point per presence pattern, sized by orthogroup count.\nFilled points = all-or-nothing candidates. None survive FDR — candidates, not biomarkers.",
                          n_holo, n_mono),
       x = "Presence-fraction difference  (holo - mono)",
       y = "-log10(Fisher p)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "pa_volcano.png"), p, width = 9, height = 5.5, dpi = 200)
ggsave(file.path(out_dir, "pa_volcano.pdf"), p, width = 9, height = 5.5)
cat("\ndone -> results/comparative/holocentricity/{pa_screen.tsv, pa_volcano.{png,pdf}}\n")
