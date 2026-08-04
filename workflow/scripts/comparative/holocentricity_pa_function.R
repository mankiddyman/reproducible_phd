#!/usr/bin/env Rscript
# =====================================================================
# Function transfer for the STRICT P/A holocentricity candidates.
# Holo-gain corner: present in all 3 holo, 0 mono  (a==3 & b==0; 105 OGs)
# Holo-loss corner: present in all 5 mono, 0 holo  (a==0 & b==5;  13 OGs)
# For each corner, pull member genes -> FANTASIA GO -> tally which
# functions appear (n OGs carrying term, n genes, n species).
# Join key: FANTASIA query_accession has a .tN transcript suffix; strip
# it to match Orthogroups.tsv gene tokens (verified ~100% overlap).
# This is raw function CONTENT, not enrichment-vs-background (next step).
# =====================================================================
suppressWarnings(suppressPackageStartupMessages({
  library(data.table); library(stringr)
}))
repo <- normalizePath(".")
of_dir <- file.path(repo, "results/comparative/orthofinder/out/Results_drosera/Orthogroups")
out_dir <- file.path(repo, "results/comparative/holocentricity")

# ---- 1) candidate OGs from the P/A screen ----
pa <- fread(file.path(out_dir, "pa_screen.tsv"))
cand <- pa[(a == 3 & b == 0) | (a == 0 & b == 5),
           .(Orthogroup, a, b,
             corner = fifelse(a == 3 & b == 0, "holo-gain", "holo-loss"))]
cat("Candidates:", nrow(cand),
    "(holo-gain", cand[corner=="holo-gain", .N], ", holo-loss", cand[corner=="holo-loss", .N], ")\n\n")

# ---- 2) melt candidate OG -> (species, gene) from Orthogroups.tsv ----
og <- fread(file.path(of_dir, "Orthogroups.tsv"))
sp_cols <- setdiff(names(og), "Orthogroup")
ogc <- og[Orthogroup %in% cand$Orthogroup]
long <- melt(ogc, id.vars = "Orthogroup", measure.vars = sp_cols,
             variable.name = "species", value.name = "genes")[genes != ""]
# split comma-separated gene lists, trim whitespace
long <- long[, .(gene = trimws(unlist(strsplit(genes, ",")))), by = .(Orthogroup, species)]
long <- merge(long, cand[, .(Orthogroup, corner)], by = "Orthogroup")
cat("Candidate member genes:", nrow(long), "across", uniqueN(long$species), "species\n\n")

# ---- 3) load FANTASIA GO for all species, strip .tN ----
fan_files <- file.path(repo, "results", sp_cols, "annotation", "function", "fantasia",
                       paste0(sp_cols, ".results.csv"))
fan <- rbindlist(lapply(seq_along(sp_cols), function(i) {
  f <- fan_files[i]; if (!file.exists(f)) return(NULL)
  d <- fread(f)
  d[, species := sp_cols[i]]
  d[, gene := sub("\\.t[0-9]+$", "", query_accession)]
  d[, .(species, gene, go_id, category, go_description, reliability_index)]
}), use.names = TRUE)
fan <- unique(fan[!is.na(go_id) & go_id != ""])
cat("FANTASIA GO rows loaded:", nrow(fan), "\n\n")

# ---- 4) join candidate genes -> GO ----
hit <- merge(long, fan, by = c("species", "gene"))
cat("Candidate gene-GO rows:", nrow(hit), "\n")
cat("Candidate genes WITH any GO:", uniqueN(hit$gene), "/", uniqueN(long$gene), "\n\n")
fwrite(hit, file.path(out_dir, "pa_strict_gene_go.tsv"), sep = "\t")

# ---- 5) per-corner GO tally ----
tally_corner <- function(cn) {
  h <- hit[corner == cn]
  if (!nrow(h)) return(data.table())
  h[, .(n_OGs = uniqueN(Orthogroup),
        n_genes = uniqueN(paste(species, gene)),
        n_species = uniqueN(species)),
    by = .(go_id, category, go_description)][order(-n_OGs, -n_genes)]
}
for (cn in c("holo-gain", "holo-loss")) {
  tb <- tally_corner(cn)
  fwrite(tb, file.path(out_dir, paste0("pa_strict_go_", sub("-", "_", cn), ".tsv")), sep = "\t")
  cat("=================================================================\n")
  cat(sprintf("%s corner: %d candidate OGs, %d GO terms\n", cn,
              cand[corner==cn, .N], nrow(tb)))
  cat("Top 25 GO terms by number of candidate OGs carrying them:\n")
  cat("=================================================================\n")
  print(head(tb[, .(go_id, cat = category,
                    description = str_trunc(go_description, 50),
                    n_OGs, n_genes, n_species)], 25))
  # also show how many candidate OGs have NO functional annotation at all
  annotated_ogs <- uniqueN(hit[corner == cn]$Orthogroup)
  cat(sprintf("\n(%d of %d %s OGs carry >=1 GO term; %d unannotated)\n\n",
              annotated_ogs, cand[corner==cn, .N], cn,
              cand[corner==cn, .N] - annotated_ogs))
}
cat("done -> pa_strict_{gene_go, go_holo_gain, go_holo_loss}.tsv\n")


# ---- 6) inspect candidate OGs with chromosome/centromere/genome-biology GO ----
# keyword grep on go_description — broad net for anything chromatin/spindle/
# centromere/division/genome-maintenance related, for manual eyeballing.
chrom_kw <- paste(c(
  "chromosome", "chromatin", "centromere", "kinetochore", "spindle",
  "microtubule", "cohesin", "condensin", "meiosis", "meiotic", "mitosis",
  "mitotic", "cell division", "cell cycle", "nucleosome", "histone",
  "telomere", "DNA replication", "DNA repair", "recombination",
  "sister chromatid", "segregation", "centrosome", "tubulin"
), collapse = "|")

chrom_hits <- hit[grepl(chrom_kw, go_description, ignore.case = TRUE)]

cat("\n=================================================================\n")
cat("Candidate OGs with chromosome/centromere/genome-biology GO terms\n")
cat("=================================================================\n")
if (!nrow(chrom_hits)) {
  cat("none found\n")
} else {
  # one row per (corner, OG, GO) — collapse genes so it's readable
  summ <- unique(chrom_hits[, .(corner, Orthogroup, go_id, category,
                                go_description)])[order(corner, Orthogroup)]
  print(summ, nrow = Inf)
  cat(sprintf("\n%d (OG, GO) pairs across %d OGs (%d holo-gain, %d holo-loss)\n",
              nrow(summ), uniqueN(summ$Orthogroup),
              uniqueN(summ[corner=="holo-gain", Orthogroup]),
              uniqueN(summ[corner=="holo-loss", Orthogroup])))
  fwrite(summ, file.path(out_dir, "pa_strict_chrom_terms.tsv"), sep = "\t")
}

