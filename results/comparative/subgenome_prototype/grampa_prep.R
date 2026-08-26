#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# GRAMPA input. Tip labels MUST be [unique gene ID]_[species label] and the tree
# must be ROOTED and BIFURCATING (polytomies are filtered by GRAMPA).
# Species labels have no underscores so the split is unambiguous.
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ape) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
dir.create("grampa", showWarnings = FALSE)

meta <- read_tsv("wgd/tip_meta.tsv", show_col_types = FALSE)
key <- function(x) gsub("@", "_", gsub("['\"]", "", x))
sp  <- c(Nepenthes_gracilis = "NEP", Dionaea_muscipula = "DIO",
         Drosera_regia = "REG", Drosera_binata = "BIN")

out <- character(0)
n_tot <- n_noroot <- n_poly <- 0
for (f in list.files("wgd/tre", "\\.tre$", full.names = TRUE)) {
  n_tot <- n_tot + 1
  tr <- tryCatch(read.tree(f), error = function(e) NULL)
  if (is.null(tr) || Ntip(tr) < 4) next
  m <- meta[match(key(tr$tip.label), key(meta$tip)), ]
  if (any(is.na(m$genome))) next
  nep <- tr$tip.label[m$genome == "Nepenthes_gracilis"]
  if (length(nep) != 1) { n_noroot <- n_noroot + 1; next }
  tr <- tryCatch(root(tr, nep, resolve.root = TRUE), error = function(e) NULL)
  if (is.null(tr)) { n_noroot <- n_noroot + 1; next }
  if (!is.binary(tr)) { n_poly <- n_poly + 1; next }
  m <- meta[match(key(tr$tip.label), key(meta$tip)), ]
  # gene ID: strip species prefix, keep unique gene name; species: short code
  gid <- sub("^[^@]*@", "", gsub("['\"]", "", m$tip))
  gid <- gsub("[^A-Za-z0-9.]", "-", gid)                  # no stray underscores
  lab <- paste0(gid, "_", unname(sp[m$genome]))
  if (anyDuplicated(lab)) next                            # GRAMPA needs unique tips
  tr$tip.label <- lab
  tr$node.label <- NULL; tr$edge.length <- NULL
  out <- c(out, write.tree(tr))
}
writeLines(out, "grampa/genetrees.txt")
writeLines("(NEP,(DIO,(REG,BIN)));", "grampa/species_tree.txt")
writeLines("(NEP,(REG,(DIO,BIN)));", "grampa/species_tree_alt.txt")

cat(sprintf("scanned %d | wrote %d | dropped: no single NEP %d, polytomy %d\n",
            n_tot, length(out), n_noroot, n_poly))
cat("example tip labels:\n")
print(head(strsplit(out[1], "[(),:]")[[1]][nzchar(strsplit(out[1], "[(),:]")[[1]])], 6))
