#!/usr/bin/env Rscript
# tpxl6_miniprot_summary.R -- reads the miniprot alignments written by tpxl6_miniprot.sh and calls each TPXL6 locus
#
# For every genome finished so far (a {species}.done file):
#   - every alignment of the two D. regia queries: TPXL6 (regia_g4881.t1) and its closest paralog from the
#     TPXL1/TPXL5 group (regia_g4900.t1)
#   - a locus hit by the TPXL6 query counts as TPXL6 only if the paralog query does not align there with a
#     higher score; otherwise it is labelled "paralog (TPXL1/5 group)"
#   - each TPXL6 locus: identity, query coverage, frameshifts, in-frame stops, exon count, intron lengths
#       intact     = >= 95% of the query aligned, no frameshift, no in-frame stop
#       disrupted  = any frameshift or in-frame stop
#       partial    = < 95% aligned, no frameshift or stop
#   - audit C1 (is the D. binata copy a pseudogene?) answered from the binata loci
# Writes tpxl6_miniprot_loci.tsv (all alignments) next to the GFFs.
# usage: Rscript tpxl6_miniprot_summary.R MINIPROT_DIR SPECIES_CSV
suppressPackageStartupMessages(library(data.table)); options(width = 200)
a <- commandArgs(TRUE)
D <- if (length(a) >= 1) a[1] else "results/comparative/holocentricity/tpx2/miniprot"
SPCSV <- if (length(a) >= 2) a[2] else "config/species.csv"
Q6 <- "regia_g4881.t1"; Q5 <- "regia_g4900.t1"
hr <- function(t) cat("\n==", t, strrep("=", max(3, 92 - nchar(t))), "\n")

ql <- local({ x <- readLines(file.path(D, "queries.faa")); h <- grepl("^>", x)
  setNames(tapply(nchar(x[!h]), cumsum(h)[!h], sum), sub("^>", "", x[h])) })
cen <- if (file.exists(SPCSV)) { s <- fread(SPCSV, select = c("species_id", "centromere")); setNames(s$centromere, s$species_id) } else character()
gen <- fread(file.path(D, "genomes.tsv"), header = FALSE, col.names = c("species", "genome"))
done <- gen[file.exists(file.path(D, paste0(species, ".done")))]
cat("genomes finished:", nrow(done), "of", nrow(gen), "\n")
if (nrow(gen) > nrow(done)) cat("  not yet:", paste(setdiff(gen$species, done$species), collapse = ", "), "\n")
if (!nrow(done)) quit(save = "no")

attr1 <- function(x, k, num = TRUE) { v <- ifelse(grepl(paste0("(^|;)", k, "="), x), sub(paste0(".*(^|;)", k, "=([^;]*).*"), "\\2", x), NA)
                                      if (num) as.numeric(v) else v }
rd <- function(sp) {
  x <- readLines(file.path(D, paste0(sp, ".gff"))); x <- x[!startsWith(x, "#") & nzchar(x)]
  if (!length(x)) return(NULL)                                   # no alignment at all in this genome
  g <- fread(text = x, header = FALSE, sep = "\t", quote = "",
             col.names = c("contig", "src", "type", "start", "end", "score", "strand", "phase", "attr"))
  m <- g[type == "mRNA"]; cds <- g[type == "CDS"]
  m[, `:=`(id = attr1(attr, "ID", FALSE), identity = attr1(attr, "Identity"), frameshift = attr1(attr, "Frameshift"),
           stop = attr1(attr, "StopCodon"), target = attr1(attr, "Target", FALSE))]
  m[is.na(frameshift), frameshift := 0]; m[is.na(stop), stop := 0]
  m[, c("query", "qs", "qe") := tstrsplit(target, " ")]
  m[, `:=`(qs = as.integer(qs), qe = as.integer(qe))]
  m[, qcov := round((qe - qs + 1) / ql[query], 3)]
  cds[, parent := attr1(attr, "Parent", FALSE)]
  ex <- cds[order(parent, start), .(n_exons = .N,
                                    introns_bp = if (.N > 1) paste(start[-1] - end[-.N] - 1, collapse = ",") else "-"), by = parent]
  m <- merge(m, ex, by.x = "id", by.y = "parent", all.x = TRUE)
  m[, species := sp]
  m[, .(species, query, contig, start, end, strand, span_kb = round((end - start + 1) / 1000, 1), score = as.numeric(score),
        identity, qcov, frameshift, stop, n_exons, introns_bp)]
}
H <- rbindlist(lapply(done$species, rd))
if (!nrow(H)) { cat("no alignments at all\n"); quit(save = "no") }
H[, centromere := fifelse(species %in% names(cen), cen[species], "outgroup")]
# assign TPXL6-query loci: paralog if the TPXL1/5-group query scores higher on an overlapping interval
H[, assigned := NA_character_]
for (i in which(H$query == Q6)) {
  o <- H[query == Q5 & species == H$species[i] & contig == H$contig[i] & start <= H$end[i] & end >= H$start[i]]
  H$assigned[i] <- if (nrow(o) && max(o$score) > H$score[i]) "paralog (TPXL1/5 group)" else "TPXL6"
}
H[query == Q5, assigned := "TPXL1/5-group query"]
H[, call := fifelse(assigned != "TPXL6", "",
            fifelse(frameshift > 0 | stop > 0, "disrupted", fifelse(qcov >= 0.95, "intact", "partial")))]
setorder(H, species, query, -score)
fwrite(H, file.path(D, "tpxl6_miniprot_loci.tsv"), sep = "\t")

hr("TPXL6 loci per genome (query: D. regia TPXL6, 139 aa)")
show <- c("contig", "start", "span_kb", "identity", "qcov", "frameshift", "stop", "n_exons", "introns_bp", "call")
for (sp in done$species) {
  x <- H[species == sp & query == Q6 & assigned == "TPXL6"]
  cat(sprintf("\n%s (%s): %d TPXL6 locus/loci\n", sp, if (sp %in% names(cen)) cen[[sp]] else "outgroup", nrow(x)))
  if (nrow(x)) print(x[, ..show], row.names = FALSE)
  p <- H[species == sp & query == Q6 & assigned != "TPXL6"]
  if (nrow(p)) cat(sprintf("   + %d TPXL6-query hit(s) where the paralog query fits better (best identity %.2f, coverage %.2f)\n",
                           nrow(p), max(p$identity), p$qcov[which.max(p$score)]))
  if (!nrow(x) && !nrow(p)) cat("   no alignment of the TPXL6 query at all\n")
}
hr("audit C1: D. binata")
b <- H[species == "Drosera_binata" & query == Q6 & assigned == "TPXL6"]
if (!"Drosera_binata" %in% done$species) cat("binata not finished yet\n") else if (!nrow(b)) cat("no TPXL6 locus in binata\n") else {
  cat(sprintf("%d binata TPXL6 locus/loci: %s\n", nrow(b), paste(sprintf("%s %s (%s)", b$contig, b$call, b$introns_bp), collapse = "; ")))
  cat(if (all(b$call == "intact")) "=> intact reading frame on every copy: no frameshift, no in-frame stop, >= 95% of the regia protein aligned.\n   The 'pseudogene' calls (tblastn frames) do not hold.\n"
      else "=> at least one copy is disrupted or partial: see the table above before choosing a figure.\n")
}
