#!/usr/bin/env Rscript
# =============================================================================
# cnv_ploidy.R  --  Copy-number comparison of holocentric vs monocentric Drosera
#
# Question: among gene families (orthogroups) present across all species, do
# holocentric species carry consistently DIFFERENT copy numbers than
# monocentric species?
#
# Design notes / caveats (read before trusting output):
#   - This is a COMPARATIVE copy-number screen on OrthoFinder gene counts.
#   - 3 holo vs 5 mono. No per-OG p-values (n is small, group structure is
#     confounded by phylogeny); candidates are called by STRICT SEPARATION
#     + effect size, NOT significance. Treat as hypothesis-generating.
#   - Raw gene counts conflate biology with (a) genome ploidy and (b)
#     annotation completeness. We correct ploidy explicitly (below). We do
#     NOT correct annotation depth -- flagged as a residual confound.
# =============================================================================

suppressWarnings(suppressPackageStartupMessages({
  library(data.table)
}))

# ---- paths ------------------------------------------------------------------
repo    <- normalizePath(".")
gc_file <- file.path(repo,
  "results/comparative/orthofinder/out/Results_drosera/Orthogroups/Orthogroups.GeneCount.tsv")
out_dir <- file.path(repo, "results/comparative/holocentricity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- species + centromere type ---------------------------------------------
# EDIT here to add/remove species (e.g. scorpioides when it lands).
cent <- c(
  Drosera_paradoxa   = "holo",
  Drosera_regia      = "holo",
  Drosera_roseana    = "holo",
  Drosera_aliciae    = "mono",
  Drosera_binata     = "mono",
  Drosera_capensis   = "mono",
  Drosera_filiformis = "mono",
  Drosera_tokaiensis = "mono"
)
species <- names(cent)
holo_sp <- names(cent)[cent == "holo"]
mono_sp <- names(cent)[cent == "mono"]
cat("HOLO:", paste(holo_sp, collapse = ", "), "\n")
cat("MONO:", paste(mono_sp, collapse = ", "), "\n\n")

# ---- load gene counts, restrict to our species ------------------------------
gc <- fread(gc_file)
if ("Total" %in% names(gc)) gc[, Total := NULL]
stopifnot(all(species %in% names(gc)))
gc <- gc[, c("Orthogroup", species), with = FALSE]

# ---- keep only orthogroups present (>0) in ALL species ----------------------
# CNV is about copy number among UNIVERSAL families. An OG absent in some
# species is a presence/absence event (different analysis), not copy number.
present_all <- gc[, rowSums(.SD > 0) == length(species), .SDcols = species]
cnv <- gc[present_all == TRUE]
cat("OGs total:", nrow(gc), "| present in all", length(species), "species:",
    nrow(cnv), "\n\n")

# =============================================================================
# PLOIDY NORMALISATION  (the core idea)
# -----------------------------------------------------------------------------
# PROBLEM: raw copy counts are not comparable across species because each
# genome sits at a different baseline ploidy / assembly redundancy:
#   - a PHASED diploid assembly reports ~2 copies of every single-copy gene
#     (both haplotypes kept separate)              -> baseline 2
#   - a COLLAPSED / diploidised genome reports ~1  -> baseline 1
#   - a recent POLYPLOID (e.g. unreduced PUTG assembly) reports ~3+
#                                                   -> baseline 3
#   If we don't correct this, a phased-diploid holo will look "expanded"
#   relative to a collapsed mono everywhere, purely as an artefact.
#
# SOLUTION: infer each species' baseline EMPIRICALLY from its own data.
#   Most gene families are single-copy-per-genome-complement, so the MOST
#   COMMON copy-number value (the MODE of the distribution across all
#   universal OGs) reflects that species' baseline ploidy. We read the
#   baseline off the data instead of assuming it.
#
#   We then divide every species' counts by ITS OWN mode, so after
#   normalisation "1.0" means "at this species' baseline" and values >1 /
#   <1 are genuine expansion / contraction relative to that baseline.
#
# WHY THIS IS DEFENSIBLE: the inferred modes recapitulate ploidy we already
#   know from the assemblies (paradoxa=2 phased, regia=1 collapsed,
#   binata=2 phased, ...), and flag known polyploids (PUTG assemblies =3).
#   So the method is validated against ground truth where we have it.
#
# RESIDUAL CAVEAT: a clean mode of 1 can mean "truly collapsed/diploidised"
#   OR "under-annotated" -- the mode cannot distinguish these. And species
#   with a smeared / bimodal distribution (low modal fraction) have an
#   unreliable baseline (e.g. capensis ~bimodal from incomplete WGD
#   diploidisation). modal_frac is reported so you can judge per species.
# =============================================================================

# infer ploidy = mode of the copy-number distribution, per species
infer_ploidy <- function(v) {
  tb <- table(factor(v, levels = 1:12))
  as.integer(names(which.max(tb)))
}
ploidy <- sapply(species, function(s) infer_ploidy(cnv[[s]]))

# modal_frac = fraction of OGs sitting at the mode; high = clean baseline,
# low (<~0.4) = smeared/bimodal -> baseline unreliable for that species.
modal_frac <- sapply(species, function(s) {
  v <- cnv[[s]]; tb <- table(factor(v, levels = 1:12))
  round(max(tb) / sum(tb), 2)
})

cat("=== inferred baseline ploidy (mode of copy-number distribution) ===\n")
print(data.table(species = species,
                 ploidy_mode = ploidy[species],
                 modal_frac  = modal_frac[species],
                 type = cent[species]))
cat("\n# modal_frac > ~0.5 = clean single peak (baseline trustworthy)\n")
cat("# modal_frac < ~0.4 = smeared/bimodal (baseline shaky -> caveat that species)\n\n")

# OPTIONAL manual override of a baseline (e.g. capensis bimodal between 1 and 2):
# ploidy["Drosera_capensis"] <- 1.5

# apply: divide each species' counts by its own baseline -> 1.0 == baseline
norm <- copy(cnv)
for (s in species) norm[[s]] <- cnv[[s]] / ploidy[s]

# sanity: per-species mean of normalised counts. Won't be exactly 1 (right
# tails of real expansions pull the mean up), but should be in the same
# ballpark across species -- big differences flag a residual depth confound.
cat("=== per-species mean copy number after ploidy-normalisation ===\n")
print(round(sapply(species, function(s) mean(norm[[s]])), 2))
cat("\n")

# =============================================================================
# CANDIDATE CALLING  (strict separation + effect size, no overlap allowed)
# =============================================================================
# group summaries on NORMALISED counts
norm[, holo_mean := rowMeans(.SD), .SDcols = holo_sp]
norm[, mono_mean := rowMeans(.SD), .SDcols = mono_sp]
norm[, holo_min  := do.call(pmin, .SD), .SDcols = holo_sp]
norm[, holo_max  := do.call(pmax, .SD), .SDcols = holo_sp]
norm[, mono_min  := do.call(pmin, .SD), .SDcols = mono_sp]
norm[, mono_max  := do.call(pmax, .SD), .SDcols = mono_sp]

# STRICT SEPARATION: every holo strictly beats every mono (or vice versa).
# This is a hard, single-species-proof bar: one noisy species breaks the call.
# High specificity, lower sensitivity -- the right trade-off for a talk.
norm[, holo_expanded   := holo_min > mono_max]
norm[, holo_contracted := holo_max < mono_min]
norm[, direction := fcase(holo_expanded,   "holo-expanded",
                          holo_contracted, "holo-contracted",
                          default = "mixed")]

# RANKING STATISTIC: separation margin = log2(weakest-in-up-group /
#   strongest-in-down-group). Using min/max (not means) means a single
#   species with a huge count CANNOT inflate the rank -- an OG only scores
#   high if the WHOLE group clears the other group with room to spare.
norm[, margin_log2 := fifelse(holo_expanded,
                              log2((holo_min + 1e-6) / (mono_max + 1e-6)),
                       fifelse(holo_contracted,
                              log2((mono_min + 1e-6) / (holo_max + 1e-6)),
                              NA_real_))]

# mean-based fold change, kept for reference (NOT used for ranking -- it is
# the statistic that lets one species dominate, which we are avoiding)
norm[, log2fc_mean := log2((holo_mean + 1e-6) / (mono_mean + 1e-6))]

# EVENNESS: how single-species-driven is an OG within the expanded group?
#   ratio of max to min member in the "up" group. ~1 = all members similar
#   (trustworthy); high = one species dominates (treat with suspicion).
norm[, holo_evenness := round(holo_max / (holo_min + 1e-6), 1)]

# ---- report ----------------------------------------------------------------
cat("=== strict separation candidates (ploidy-normalised) ===\n")
cat("holo-expanded   (all holo > all mono):", norm[holo_expanded   == TRUE, .N], "\n")
cat("holo-contracted (all holo < all mono):", norm[holo_contracted == TRUE, .N], "\n")
cat("mixed:", norm[direction == "mixed", .N], "\n\n")

# save full table
fwrite(norm, file.path(out_dir, "cnv_ploidy_screen.tsv"), sep = "\t")

# pretty per-species columns for printing
show_cols <- function(dt, n = 30) {
  dt[, c("Orthogroup", species,
         "margin_log2", "holo_evenness", "log2fc_mean"), with = FALSE
     ][, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)][seq_len(min(n, .N))]
}

cat("=== top holo-expanded by separation margin (single-species-proof) ===\n")
print(show_cols(norm[holo_expanded == TRUE][order(-margin_log2)], 30))

cat("\n=== holo-contracted by separation margin ===\n")
print(show_cols(norm[holo_contracted == TRUE][order(-margin_log2)], 30))

cat("\nWrote: ", file.path(out_dir, "cnv_ploidy_screen.tsv"), "\n")


suppressWarnings(suppressPackageStartupMessages(library(data.table)))
repo    <- normalizePath(".")
out_dir <- file.path(repo, "results/comparative/holocentricity")
of_dir  <- file.path(repo, "results/comparative/orthofinder/out/Results_drosera/Orthogroups")

# ---- rebuild OG -> GO map (join key = species+gene, never gene alone) ----
og_tsv <- fread(file.path(of_dir, "Orthogroups.tsv"))
sp_all <- setdiff(names(og_tsv), "Orthogroup")
long <- melt(og_tsv, id.vars="Orthogroup", measure.vars=sp_all,
             variable.name="species", value.name="genes")[genes != ""]
long <- long[, .(gene=trimws(unlist(strsplit(genes, ",")))), by=.(Orthogroup, species)]
long[, species := as.character(species)]
fan <- rbindlist(lapply(sp_all, function(s){
  f <- file.path(repo,"results",s,"annotation","function","fantasia",paste0(s,".results.csv"))
  if(!file.exists(f)) return(NULL)
  d <- fread(f); d[, species := s]; d[, gene := sub("\\.t[0-9]+$","",query_accession)]
  d[, .(species, gene, go_id, category, go_description)]
}), use.names=TRUE)
fan <- unique(fan[!is.na(go_id) & go_id!=""])
og_go <- unique(merge(long, fan, by=c("species","gene"))[, .(Orthogroup, go_id, category, go_description)])
cat("OG->GO pairs:", nrow(og_go), "(expect ~206465)\n")

# ---- reload the screen + define OGs of interest ----
scr <- fread(file.path(out_dir, "cnv_ploidy_screen.tsv"))

# the 3 credible high-margin holo-expanded
high_margin <- c("OG0002011","OG0000667","OG0001376")

# GO terms for the 3 of interest
cat("\n=== GO terms: 3 high-margin holo-expanded OGs ===\n")
for(og in high_margin){
  terms <- og_go[Orthogroup==og, sort(unique(go_description))]
  cat("\n---", og, "(margin", round(scr[Orthogroup==og, margin_log2],2),
      "| holo", paste(round(unlist(scr[Orthogroup==og, .(Drosera_paradoxa,Drosera_regia,Drosera_roseana)]),1),collapse="/"),
      ") ---\n")
  cat(paste(terms, collapse="; "), "\n")
}
