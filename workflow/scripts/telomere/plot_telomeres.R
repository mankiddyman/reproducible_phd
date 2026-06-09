#!/usr/bin/env Rscript
# Telomere / gap / contig visualisation for a scaffolded assembly.
# Standalone helper for the manual telomere-flipping curation loop (NOT a
# snakemake rule). Invoked by teloflip.sh, but runnable directly:
#
#   module load genespace
#   Rscript plot_telomeres.R <assembly.fa> <n_scaffolds> [hold_seconds]
#
# Opens an X11 window (over ssh -X) showing telomere positions on
# scaffold_1..n so you can see which need flipping so telomeres sit at the
# terminals. Window is held until you press Enter, or for <hold_seconds>.

suppressPackageStartupMessages({
  library(GENESPACE)
  library(Biostrings)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: Rscript plot_telomeres.R <assembly.fa> <n_scaffolds> [hold_seconds]")
}
asm       <- args[[1]]
n_scaff   <- as.integer(args[[2]])
hold_secs <- if (length(args) >= 3) as.numeric(args[[3]]) else NA

if (!file.exists(asm)) stop(sprintf("assembly not found: %s", asm))
message(sprintf("Assembly:  %s", asm))
message(sprintf("Scaffolds: scaffold_1..%d", n_scaff))

dnass <- readDNAStringSet(asm)

# Subset to chromosome-scale scaffolds BEFORE telomere-finding. Running
# find_contigsGapsTelos on all (30k+) unplaced contigs hangs for many minutes;
# we only plot scaffold_1..n anyway, so restrict up front (seconds not minutes).
valid_scaffolds_pre <- paste0("scaffold_", seq_len(n_scaff))
dnass <- dnass[names(dnass) %in% valid_scaffolds_pre]
message(sprintf("Subset to %d chromosome scaffolds before telomere-finding", length(dnass)))

# Plant telomere repeat (Arabidopsis-type), both strands.
teloKmers <- c("TTTAGGG", "CCCTAAA")

result <- find_contigsGapsTelos(
  dnass            = dnass,
  teloKmers        = teloKmers,
  minContigGapSize = 100,    # min N-run size to call a gap
  maxDistBtwTelo   = 20,     # max gap between neighbouring telo kmers
  minTeloSize      = 200,    # min telomere kmer-cluster size
  minTeloDens      = 0.75,   # min telomere kmer density (0-1)
  minChrSize       = 0,      # min scaffold size to include
  maxDist2end      = 10000,  # max distance to chr end for a telomere call
  verbose          = TRUE
)

# Keep only scaffold_1..n (the chromosome-scale scaffolds).
valid_scaffolds <- paste0("scaffold_", seq_len(n_scaff))
filtered <- lapply(result, function(gr) gr[seqnames(gr) %in% valid_scaffolds])

palette <- colorRampPalette(c("blue", "green"))

# Open an explicit X11 device. Rscript is non-interactive, so unlike the REPL
# we must open the device ourselves AND keep the process alive, else the
# window closes the instant the script ends.
opened_x11 <- FALSE
tryCatch({
  x11(width = 10, height = 12, type = "cairo")
  opened_x11 <- TRUE
}, error = function(e) {
  message("x11() failed (", conditionMessage(e),
          ") - falling back to PNG next to the assembly.")
})

plot_contigs(cgt = filtered, nColors = 4, palette = palette)

# Print the >1kb telomere table - handy for deciding flips.
telo_df <- as.data.frame(result$telo)
big <- telo_df[telo_df$width > 1000, ]
if (nrow(big)) {
  message("\nTelomere clusters > 1kb:")
  print(big[order(big$seqnames, big$start), c("seqnames","start","end","width")])
} else {
  message("\n(no telomere clusters > 1kb)")
}

if (opened_x11) {
  if (is.na(hold_secs)) {
    cat("\n[Enter] to close the plot window... ")
    invisible(readLines("/dev/tty", n = 1))   # blocks under Rscript; readline() does not
  } else {
    Sys.sleep(hold_secs)
  }
} else {
  png_path <- sub("\\.(fa|fasta)$", ".telomeres.png", asm)
  png(png_path, width = 1000, height = 1200)
  plot_contigs(cgt = filtered, nColors = 4, palette = palette)
  dev.off()
  message(sprintf("Wrote plot: %s", png_path))
}
