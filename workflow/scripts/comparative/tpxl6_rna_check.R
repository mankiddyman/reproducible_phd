#!/usr/bin/env Rscript
# tpxl6_rna_check.R -- is TPXL6 expressed in leaf RNA? Reads counted at the gene, not read off annotation tags
#
# Why: "D. binata TPXL6 has no leaf RNA" so far comes from BRAKER's evidence tags (gene_support.tsv), which only
# say that no RNA hint supported a gene model there. Here the reads are counted directly.
# Data: the STAR 2-pass alignments of each species' leaf RNA made for the annotation bake-off (rule star_align):
#   results/{species}/rna_support/star_align/{Aligned.sortedByCoord.out.bam, SJ.out.tab, Log.final.out}
#   (D. binata and D. paradoxa; species without them are reported and skipped). The BAM is indexed if needed.
# Loci: every TPXL6 copy found by miniprot (tpxl6_miniprot_loci.tsv, assigned TPXL6) and, as a control in the
# same library, the best-scoring locus of the TPXL1/TPXL5-group query (TPXL6's closest relatives).
#   reads      primary alignments overlapping the exons (samtools view -c -F 0x904 -M -L exons.bed). In a phased
#              assembly a read from the gene lands on one haplotype or the other, so the copies are also summed.
#   unique     the same, MAPQ 255 only (STAR's uniquely mapped reads)
#   junctions  STAR splice junctions (SJ.out.tab) within 3 bp of the copy's introns, and the reads spanning them
#   RPKM       reads per kb of exon per million mapped reads (unique + multi, from Log.final.out)
# usage: Rscript tpxl6_rna_check.R REPO MINIPROT_DIR SAMTOOLS OUT_TSV
suppressPackageStartupMessages(library(data.table)); options(width = 200)
a <- commandArgs(TRUE)
if (length(a) < 4) stop("usage: tpxl6_rna_check.R REPO MINIPROT_DIR SAMTOOLS OUT_TSV")
REPO <- a[1]; D <- a[2]; SAM <- a[3]; OUT <- a[4]
Q6 <- "regia_g4881.t1"; Q5 <- "regia_g4900.t1"
hr <- function(t) cat("\n==", t, strrep("=", max(3, 92 - nchar(t))), "\n")
attr1 <- function(x, k) ifelse(grepl(paste0("(^|;)", k, "="), x), sub(paste0(".*(^|;)", k, "=([^;]*).*"), "\\2", x), NA)
gff_of <- function(sp) {
  x <- readLines(file.path(D, paste0(sp, ".gff"))); x <- x[!startsWith(x, "#") & nzchar(x)]
  g <- fread(text = x, header = FALSE, sep = "\t", quote = "",
             col.names = c("contig", "src", "type", "start", "end", "score", "strand", "phase", "attr"))
  g[, `:=`(id = attr1(attr, "ID"), parent = attr1(attr, "Parent"), target = attr1(attr, "Target"))]
  g
}
count <- function(args) { x <- suppressWarnings(system2(SAM, args, stdout = TRUE, stderr = FALSE)); as.integer(x[length(x)]) }
L <- fread(file.path(D, "tpxl6_miniprot_loci.tsv"))
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

hr("leaf RNA at TPXL6 and at a TPXL1/5-group relative")
res <- list()
for (sp in unique(L$species)) {
  sd <- file.path(REPO, "results", sp, "rna_support", "star_align")
  bam <- file.path(sd, "Aligned.sortedByCoord.out.bam"); sjf <- file.path(sd, "SJ.out.tab"); lg <- file.path(sd, "Log.final.out")
  if (!file.exists(bam)) { cat(sprintf("%-20s no STAR alignment, skipped\n", sp)); next }
  if (!any(file.exists(paste0(bam, c(".bai", ".csi"))))) {
    cat(sprintf("%-20s indexing %s (once)\n", sp, bam)); system2(SAM, c("index", "-@", "8", bam))
  }
  lib <- NA_real_
  if (file.exists(lg)) { x <- readLines(lg)
    gv <- function(k) as.numeric(trimws(sub(".*\\|", "", grep(k, x, value = TRUE, fixed = TRUE)[1])))
    lib <- gv("Uniquely mapped reads number") + gv("Number of reads mapped to multiple loci") }
  SJ <- if (file.exists(sjf) && file.size(sjf) > 0)
          fread(sjf, header = FALSE, col.names = c("chr", "s", "e", "strand", "motif", "annot", "uniq", "multi", "ovh")) else NULL
  if (is.null(SJ)) cat(sprintf("%-20s no splice-junction file (or it is empty): junctions counted as 0\n", sp))
  sq <- sub(".*\tSN:([^\t]+).*", "\\1", grep("^@SQ", system2(SAM, c("view", "-H", bam), stdout = TRUE), value = TRUE))
  g <- gff_of(sp)
  loci <- rbind(copy(L[species == sp & query == Q6 & assigned == "TPXL6"])[, class := "TPXL6"],
                copy(L[species == sp & query == Q5][order(-score)][1])[, class := "TPXL1/5-group control"], fill = TRUE)
  for (i in seq_len(nrow(loci))) {
    r <- loci[i]; if (is.na(r$contig)) next
    if (!r$contig %in% sq) { cat(sprintf("%-20s %s is not in this BAM: alignment made on another assembly? skipped\n", sp, r$contig)); next }
    m <- g[type == "mRNA" & contig == r$contig & start == r$start & end == r$end & startsWith(target, paste0(r$query, " "))][1]
    ex <- g[type == "CDS" & parent == m$id][order(start)]
    bed <- tempfile(fileext = ".bed"); fwrite(ex[, .(contig, start - 1L, end)], bed, sep = "\t", col.names = FALSE)
    n_all <- count(c("view", "-c", "-F", "0x904", "-M", "-L", bed, bam))
    n_uni <- count(c("view", "-c", "-F", "0x904", "-q", "255", "-M", "-L", bed, bam))
    unlink(bed)
    intr <- if (nrow(ex) > 1) data.table(s = ex$end[-nrow(ex)] + 1L, e = ex$start[-1] - 1L) else data.table(s = integer(), e = integer())
    jn <- if (!is.null(SJ) && nrow(intr)) sapply(seq_len(nrow(intr)), function(k) {
      h <- SJ[chr == r$contig & abs(s - intr$s[k]) <= 3 & abs(e - intr$e[k]) <= 3]; sum(h$uniq + h$multi) }) else integer()
    exbp <- sum(ex$end - ex$start + 1L)
    res[[length(res) + 1]] <- data.table(species = sp, class = r$class, contig = r$contig, start = r$start, exons = nrow(ex),
      exon_bp = exbp, reads = n_all, unique_reads = n_uni, introns = nrow(intr), introns_with_junction = sum(jn > 0),
      junction_reads = sum(jn), library_reads = lib, RPKM = round(n_all / (exbp / 1000) / (lib / 1e6), 2))
  }
}
if (!length(res)) { cat("no species with a STAR alignment: nothing to count\n"); quit(save = "no") }
R <- rbindlist(res)
fwrite(R, OUT, sep = "\t")
print(R[, .(species, class, contig, exons, reads, unique_reads, junctions = paste0(introns_with_junction, "/", introns), junction_reads, RPKM)], row.names = FALSE)

hr("per species (copies summed)")
S <- R[, .(copies = .N, reads = sum(reads), junction_reads = sum(junction_reads),
           RPKM = round(sum(reads) / (mean(exon_bp) / 1000) / (library_reads[1] / 1e6), 2)), by = .(species, class)]
print(S, row.names = FALSE)
for (sp in unique(S$species)) {
  t6 <- S[species == sp & class == "TPXL6"]; ct <- S[species == sp & class != "TPXL6"]
  if (!nrow(t6)) next
  cat(sprintf("\n%s: TPXL6 %d reads, %d junction reads, RPKM %s", sp, t6$reads, t6$junction_reads, t6$RPKM))
  if (nrow(ct)) cat(sprintf(" | TPXL1/5-group relative %d reads, RPKM %s | TPXL6 at %.1f%% of the relative", ct$reads, ct$RPKM,
                            if (ct$RPKM > 0) 100 * t6$RPKM / ct$RPKM else NA))
}
cat("\n\nwrote", OUT, "\n")
