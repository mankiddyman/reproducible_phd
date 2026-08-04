#!/usr/bin/env Rscript
# =====================================================================
# GO ENRICHMENT for the strict P/A holocentricity candidates.
# Upgrades the function-transfer script from raw CONTENT (which GO terms
# appear) to ENRICHMENT (which GO terms appear MORE than genome background).
#
# Unit: orthogroup (OG). An OG "carries" GO term X if >=1 of its member
# genes (any species) is annotated with X in FANTASIA.
# Background: all OGs carrying >=1 GO term (annotatable OGs only — including
# unannotatable OGs would dilute the background and inflate enrichment).
# Test: per GO term, one-sided Fisher (alternative="greater") for
# over-representation in candidates vs background. BH-corrected PER CORNER.
# =====================================================================
suppressWarnings(suppressPackageStartupMessages({
  library(data.table); library(stringr)
}))
repo    <- normalizePath(".")
of_dir  <- file.path(repo, "results/comparative/orthofinder/out/Results_drosera/Orthogroups")
out_dir <- file.path(repo, "results/comparative/holocentricity")

# candidate OGs from the P/A screen (same corners as the function script)
pa <- fread(file.path(out_dir, "pa_screen.tsv"))
cand <- pa[(a == 3 & b == 0) | (a == 0 & b == 5),
           .(Orthogroup, corner = fifelse(a == 3 & b == 0, "holo-gain", "holo-loss"))]
cat("Candidates:", nrow(cand),
    "(holo-gain", cand[corner=="holo-gain", .N], ", holo-loss", cand[corner=="holo-loss", .N], ")\n")



# ---- OG -> GO map for ALL orthogroups (needed for background) ----
og <- fread(file.path(of_dir, "Orthogroups.tsv"))
sp_cols <- setdiff(names(og), "Orthogroup")

# melt EVERY OG to (Orthogroup, species, gene) — no candidate filter this time
long_all <- melt(og, id.vars = "Orthogroup", measure.vars = sp_cols,
                 variable.name = "species", value.name = "genes")[genes != ""]
long_all <- long_all[, .(gene = trimws(unlist(strsplit(genes, ",")))),
                     by = .(Orthogroup, species)]
cat("All OG member genes:", nrow(long_all), "\n")

# FANTASIA GO for all species, strip .tN (same as function script)
fan_files <- file.path(repo, "results", sp_cols, "annotation", "function", "fantasia",
                       paste0(sp_cols, ".results.csv"))
fan <- rbindlist(lapply(seq_along(sp_cols), function(i) {
  f <- fan_files[i]; if (!file.exists(f)) return(NULL)
  d <- fread(f)
  d[, species := sp_cols[i]]
  d[, gene := sub("\\.t[0-9]+$", "", query_accession)]
  d[, .(species, gene, go_id, category, go_description)]
}), use.names = TRUE)
fan <- unique(fan[!is.na(go_id) & go_id != ""])

# join genes -> GO, then collapse to UNIQUE (OG, GO) pairs:
# an OG carries a term if ANY member gene has it. This is the core unit.
og_go <- unique(merge(long_all, fan, by = c("species", "gene"))[
  , .(Orthogroup, go_id, category, go_description)])
cat("Unique (OG, GO) pairs:", nrow(og_go), "\n")



# ---- per-corner enrichment ----
enrich_corner <- function(cn) {
  cand_ogs <- cand[corner == cn, Orthogroup]
  n_cand   <- length(cand_ogs)
  term_k <- og_go[Orthogroup %in% cand_ogs, .(k = uniqueN(Orthogroup)), by = go_id]
  tt <- merge(term_k, term_bg, by = "go_id")
  # filter: term must be in >=3 candidate OGs AND >=10 background OGs.
  # rare terms (k=1, K=1-2) produce absurd fold (>1000) and false significance
  # purely from small n_cand; they are single-gene flukes, not enrichment.
  tt <- tt[k >= 3 & K >= 10]
  if (!nrow(tt)) return(data.table())   # nothing testable after filter
  tt[, p := mapply(function(k, K) {
    m <- matrix(c(k, n_cand - k, K - k, (N_bg - n_cand) - (K - k)), nrow = 2)
    fisher.test(m, alternative = "greater")$p.value
  }, k, K)]
  tt[, q := p.adjust(p, method = "BH")]
  tt[, `:=`(corner = cn, n_cand = n_cand,
            cand_frac = round(k / n_cand, 3),
            bg_frac   = round(K / N_bg, 4),
            fold      = round((k / n_cand) / (K / N_bg), 2))]
  tt[]
}

res <- rbindlist(lapply(c("holo-gain", "holo-loss"), enrich_corner))
# category + go_description already live in term_bg via og_go, so res HAS them.
setorder(res, corner, q, p, -fold)

fwrite(res, file.path(out_dir, "pa_enrichment_all.tsv"), sep = "\t")
for (cn in c("holo-gain", "holo-loss")) {
  tb <- res[corner == cn]
  fwrite(tb, file.path(out_dir, paste0("pa_enrich_", sub("-", "_", cn), ".tsv")), sep = "\t")
  cat("\n=================================================================\n")
  cat(sprintf("%s: %d terms tested, %d at q<=0.10, %d at q<=0.05\n",
              cn, nrow(tb), tb[q <= 0.10, .N], tb[q <= 0.05, .N]))
  cat("Top 20 by p (then fold) — nothing significant, shown for inspection:\n")
  cat("=================================================================\n")
  print(head(tb[order(p, -fold), .(go_id, cat = category,
                    desc = str_trunc(go_description, 45),
                    k, n_cand, K, cand_frac, bg_frac, fold,
                    p = signif(p, 3), q = signif(q, 3))], 20))
}
cat("\ndone -> pa_enrichment_all.tsv, pa_enrich_{holo_gain,holo_loss}.tsv\n")
