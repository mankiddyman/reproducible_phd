#!/usr/bin/env Rscript
# GENESPACE v1.3.1. Called by rule genespace_run (workflow/rules/synteny.smk).
# Args: <wd> <path2mcscanx> <nCores> <refGenome> <genomeID1,genomeID2,...>
# wd must already contain bed/ and peptide/; GENESPACE writes its outputs there.
suppressMessages(library(GENESPACE))

args    <- commandArgs(trailingOnly = TRUE)
wd      <- normalizePath(args[1], mustWork = TRUE)
mcscanx <- normalizePath(args[2], mustWork = TRUE)
ncores  <- as.integer(args[3])
refGen  <- args[4]
genomeIDs <- strsplit(args[5], ",")[[1]]

cat("wd:      ", wd, "\nmcscanx: ", mcscanx, "\nnCores:  ", ncores,
    "\nref:     ", refGen, "\ngenomes: ", paste(genomeIDs, collapse = ", "), "\n")

# fails loudly on a malformed bed/peptide pair rather than deep inside the run
for (g in genomeIDs) check_annotFiles(filepath = wd, genomeIDs = c(g))

gpar <- init_genespace(wd = wd, path2mcscanx = mcscanx,
                       genomeIDs = genomeIDs, nCores = ncores)
out  <- run_genespace(gpar, overwrite = TRUE)

# --- white-background riparian (supervisor-facing) ---
ggthemes  <- ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(c("#D61F27", "#F9AC60", "#FCF7BF", "#ADD8E7", "#307BB5"))

pdf(file = file.path(wd, "riparian_plot_white.pdf"), width = 8, height = 5)
plot_riparian(gsParam = out, palette = customPal, addThemes = ggthemes,
              chrFill = "lightgrey", backgroundColor = NULL,
              refGenome = refGen, genomeIDs = genomeIDs)
dev.off()

# default dark version too, for comparison
pdf(file = file.path(wd, "riparian_plot_default.pdf"), width = 8, height = 5)
plot_riparian(gsParam = out, refGenome = refGen, genomeIDs = genomeIDs)
dev.off()

cat("=== genespace complete ===\n")
