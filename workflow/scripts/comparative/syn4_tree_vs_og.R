#!/usr/bin/env Rscript
# syn4_tree_vs_og.R -- do the SYN4 gene tree and the SYN4 violin count the same thing?
#
# The violin counts the genes of one OrthoFinder orthogroup (OG) per species: one protein per gene
# (the longest isoform), haplotype copies counted separately. A gene tree can disagree with that count in
# four ways, and this script tells them apart tip by tip:
#   isoforms  several tips are transcripts (.t1 .t2 ...) of ONE gene        -> the tree overcounts
#   split OG  a tip's gene sits in ANOTHER OG, or in no OG (unassigned)    -> the violin undercounts
#   version   a tip matches no gene of the OrthoFinder input               -> the two used different annotations
#   fragment  a short sequence (< half the median length), or one labelled frag/partial/pseudo
# For a tree holding only SYN4-type sequences it also roots on the non-Droseraceae tips (Arabidopsis, Beta,
# Nepenthes) and lists, per species, which genes fall in each daughter clade of the Droseraceae node.
#
# Usage: Rscript syn4_tree_vs_og.R OUT_DIR ORTHOGROUPS_DIR REPO OG_ID TREE [TREE ...]
# Writes OUT_DIR/{tips,summary}_<tree>.tsv and OUT_DIR/topology_<tree>.nwk; prints the comparison.
suppressPackageStartupMessages({ library(ape); library(data.table) })
options(width = 250)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("usage: OUT_DIR ORTHOGROUPS_DIR REPO OG_ID TREE [TREE ...]")
OUT <- args[1]; OGD <- args[2]; REPO <- args[3]; OGID <- args[4]; TREES <- args[-(1:4)]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

DROS  <- c("regia", "paradoxa", "roseana", "scorpioides", "binata", "capensis", "aliciae", "tokaiensis", "filiformis")
HOLO  <- c("regia", "paradoxa", "roseana", "scorpioides")
OTHER <- c(arabidopsis = "Arabidopsis", atha = "Arabidopsis", ath = "Arabidopsis", athal = "Arabidopsis",
           beta = "Beta", bvul = "Beta", bvulgaris = "Beta",
           nepenthes = "Nepenthes", ngra = "Nepenthes", ngracilis = "Nepenthes",
           dionaea = "Dionaea", dmus = "Dionaea", dmuscipula = "Dionaea")
OUTGROUPS <- c("Arabidopsis", "Beta", "Nepenthes")
FULLNAME <- c(setNames(paste0("Drosera_", DROS), DROS), dionaea = "Dionaea_muscipula", nepenthes = "Nepenthes_gracilis")
FAMILY_MARKERS <- "AT5G05490|AT5G40840|AT3G59550|SYN1|SYN2|SYN3|REC8|DIF1"   # the other alpha-kleisins
hr <- function(t) cat(sprintf("\n== %s\n", t))

# ---------------------------------------------------------------- OrthoFinder: every gene -> its OG
read_ogs <- function(f, unassigned) {
  if (!file.exists(f)) return(NULL)
  x <- fread(f, sep = "\t", header = TRUE, colClasses = "character", quote = "")
  long <- melt(x, id.vars = "Orthogroup", variable.name = "col", value.name = "genes", variable.factor = FALSE)
  long <- long[!is.na(genes) & genes != ""]
  long <- long[, .(gene = trimws(unlist(strsplit(genes, ",", fixed = TRUE)))), by = .(Orthogroup, col)]
  long[, unassigned := unassigned]
  long
}
ogf <- file.path(OGD, "Orthogroups.tsv")
if (!file.exists(ogf)) stop("no Orthogroups.tsv in ", OGD)
og <- rbind(read_ogs(ogf, FALSE), read_ogs(file.path(OGD, "Orthogroups_UnassignedGenes.tsv"), TRUE))
og[, sp := tolower(sub("^Drosera_", "", col))]
og[sp == "dionaea_muscipula", sp := "dionaea"]
og[sp == "nepenthes_gracilis", sp := "nepenthes"]
og[, gene_bare := sub("_(rnasupp|norna)$", "", gene)]
by_gene <- og[, .(key = paste(sp, gene), gene, Orthogroup, unassigned)]
by_bare <- og[, .(key = paste(sp, gene_bare), gene, Orthogroup, unassigned)]
setkey(by_gene, key); setkey(by_bare, key)
run_species <- sort(unique(og$sp))

gc_file <- file.path(OGD, "Orthogroups.GeneCount.tsv")
gcount <- if (file.exists(gc_file)) fread(gc_file, sep = "\t", header = TRUE) else NULL
hr(sprintf("OrthoFinder run: %s", OGD))
cat(sprintf("   species (%d): %s\n", length(run_species), paste(run_species, collapse = " ")))
cat(sprintf("   orthogroups: %s; genes: %s (%s unassigned); Orthogroups.tsv modified %s\n",
            format(uniqueN(og[unassigned == FALSE]$Orthogroup), big.mark = ","), format(nrow(og), big.mark = ","),
            format(sum(og$unassigned), big.mark = ","), format(file.mtime(ogf), "%F %H:%M")))

og_members <- og[Orthogroup == OGID & unassigned == FALSE]
violin <- setNames(integer(0), character(0))
if (!is.null(gcount) && OGID %in% gcount$Orthogroup) {
  r <- gcount[Orthogroup == OGID]
  cols <- setdiff(names(r), c("Orthogroup", "Total"))
  violin <- setNames(as.integer(unlist(r[, ..cols])), tolower(sub("^Drosera_", "", cols)))
}
hr(sprintf("%s in this run, genes per species (the violin's diamonds, if this is the run it used)", OGID))
if (!length(violin)) cat("   !! not in this run's GeneCount table: OrthoFinder was rerun and the IDs renumbered\n")
for (s in names(violin)) cat(sprintf("   %-12s %2d  %s\n", s, violin[[s]], paste(og_members[sp == s]$gene, collapse = " ")))

# ---------------------------------------------------------------- annotation: gene -> locus
gff_cache <- list()
gff_genes <- function(s) {
  if (!is.null(gff_cache[[s]])) return(gff_cache[[s]])
  f <- file.path(REPO, "results", FULLNAME[[s]], "annotation", "final", paste0(FULLNAME[[s]], ".final.gff3"))
  res <- data.table(id = character(), id_bare = character(), seqid = character(), start = integer(), end = integer())
  if (file.exists(f)) {
    g <- tryCatch(fread(cmd = sprintf("awk -F'\\t' '$3==\"gene\"' '%s'", f), sep = "\t", header = FALSE, quote = "",
                        select = c(1, 4, 5, 9), colClasses = "character"), error = function(e) data.table())
    if (nrow(g)) {
      setnames(g, c("seqid", "start", "end", "attr"))
      g[, id := sub("^(.*;)?ID=([^;]+).*$", "\\2", attr)]
      res <- g[, .(id, id_bare = sub("_(rnasupp|norna)$", "", id), seqid, start = as.integer(start), end = as.integer(end))]
    }
  }
  gff_cache[[s]] <<- res
  res
}

# ---------------------------------------------------------------- tip parsing
norm <- function(x) gsub("[^A-Za-z0-9.-]", "_", x)
species_of <- function(lab) {
  l <- tolower(lab)
  keys <- c(DROS, names(OTHER)); keys <- keys[order(-nchar(keys))]
  for (k in keys) if (grepl(paste0("^(drosera_|d_|d\\.)?", k, "([_|.: -]|$)"), l)) return(k)
  for (k in c(DROS, "arabidopsis", "nepenthes", "dionaea", "bvulgaris", "dmuscipula", "ngracilis")) if (grepl(k, l, fixed = TRUE)) return(k)
  NA_character_
}
strip_prefix <- function(lab, k) {
  if (is.na(k)) return(lab)
  out <- sub(paste0("^(?i)(drosera_|d_|d\\.)?", k, "[_|.: -]+"), "", lab, perl = TRUE)
  if (identical(out, lab)) out <- sub(paste0("(?i)^.*", k, "[_|.: -]+"), "", lab, perl = TRUE)
  sub("^[_|]+", "", out)
}
match_gene <- function(k, core) {
  none <- list(gene = NA_character_, og = NA_character_, unassigned = FALSE)
  if (is.na(k) || !(k %in% run_species)) return(none)
  x <- core
  for (step in 1:6) {
    cands <- unique(c(x, sub("\\.t[0-9]+$", "", x), sub("\\.[0-9]+$", "", x), sub("[.-](mRNA|RA|PA)[0-9]*$", "", x)))
    for (cand in cands) {
      h <- by_gene[J(paste(k, cand)), nomatch = NULL]
      if (!nrow(h)) h <- by_bare[J(paste(k, sub("_(rnasupp|norna)$", "", cand))), nomatch = NULL]
      if (nrow(h)) return(list(gene = h$gene[1], og = h$Orthogroup[1], unassigned = h$unassigned[1]))
    }
    if (grepl("[_|][^_|]+$", x)) x <- sub("[_|][^_|]+$", "", x) else break
  }
  none
}
read_fasta <- function(f) {
  l <- readLines(f, warn = FALSE)
  h <- grep("^>", l)
  if (!length(h)) return(NULL)
  ids <- sub("^>(\\S+).*$", "\\1", l[h])
  ends <- c(h[-1] - 1, length(l))
  seqs <- mapply(function(a, b) if (b > a) paste(l[(a + 1):b], collapse = "") else "", h, ends)
  data.table(id = ids, nid = norm(ids), seq = toupper(gsub("[-*. ]", "", seqs)))
}
find_alignment <- function(tree_file, labels) {
  d <- dirname(tree_file)
  stem <- sub("\\.(treefile|contree|nwk|newick|tre|tree)$", "", basename(tree_file))
  fs <- list.files(d, pattern = "\\.(fa|faa|fasta|fas|fna|aln|afa|mafft|msa|trim|trimal|clipkit)$", full.names = TRUE)
  fs <- fs[file.size(fs) < 50e6]
  if (!length(fs)) return(NULL)
  fs <- fs[order(-startsWith(basename(fs), stem))]
  best <- NULL
  for (f in fs) {
    fa <- tryCatch(read_fasta(f), error = function(e) NULL)
    if (is.null(fa)) next
    cov <- mean(norm(labels) %in% fa$nid)
    if (is.null(best) || cov > best$cov) best <- list(file = f, fa = fa, cov = cov)
  }
  if (!is.null(best) && best$cov >= 0.5) best else NULL
}
loci_proxy <- function(seqids) {
  # copies per haploid set. Chromosome-level phased assemblies: chrN_hap1 + chrN_hap2 is one locus seen twice,
  # two copies on one haplotype's chromosome are two loci. Contig-level phased (h1tg/h2tg contigs): haplotypes
  # cannot be paired by name, so take the larger of the two haplotypes' counts. Collapsed: every copy counts.
  if (!length(seqids)) return(0L)
  if (any(grepl("_hap[12]$", seqids))) {
    base <- sub("_hap[12]$", "", seqids)
    hap  <- ifelse(grepl("_hap[12]$", seqids), sub("^.*_(hap[12])$", "\\1", seqids), "one")
    dt <- data.table(base, hap)[, .N, by = .(base, hap)]
    return(as.integer(sum(dt[, .(m = max(N)), by = base]$m)))
  }
  if (all(grepl("^h[12]tg", seqids))) return(as.integer(max(table(substr(seqids, 1, 2)))))
  length(seqids)
}

# ---------------------------------------------------------------- per tree
for (tf in TREES) {
  tr <- tryCatch(read.tree(tf), error = function(e) NULL)
  if (is.null(tr)) { hr(sprintf("tree: %s -- could not be read", tf)); next }
  if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
  tag <- gsub("[^A-Za-z0-9]+", "_", sub("\\.[^.]*$", "", basename(tf)))
  labs <- tr$tip.label
  family <- any(grepl(FAMILY_MARKERS, labs, ignore.case = TRUE))
  aln <- find_alignment(tf, labs)

  T <- data.table(tip = labs)
  T[, k := vapply(tip, species_of, "")]
  T[, group := ifelse(k %in% DROS, k, ifelse(is.na(k), "unknown", OTHER[k]))]
  T[, core := unname(mapply(strip_prefix, tip, k))]
  T[, isoform := ifelse(grepl("\\.t[0-9]+", core), sub("^.*(\\.t[0-9]+).*$", "\\1", core), "")]
  m <- lapply(seq_len(nrow(T)), function(i) match_gene(T$k[i], T$core[i]))
  T[, gene := vapply(m, function(z) z$gene, "")]
  T[, og := vapply(m, function(z) z$og, "")]
  T[, unassigned := vapply(m, function(z) isTRUE(z$unassigned), TRUE)]
  T[is.na(gene) & group %in% c(OUTGROUPS, "Dionaea"), gene := sub("\\.[0-9]+$", "", sub("\\.t[0-9]+$", "", core))]
  T[, in_og := !is.na(og) & og == OGID]
  T[, len := NA_integer_][, seqgrp := NA_integer_]
  if (!is.null(aln)) {
    sq <- aln$fa$seq[match(norm(T$tip), aln$fa$nid)]
    T[, len := nchar(sq)][is.na(sq), len := NA_integer_]
    T[, seqgrp := match(sq, unique(sq[!is.na(sq)]))]
  }
  medlen <- suppressWarnings(median(T[group %in% c(DROS, "Dionaea") & !is.na(len)]$len))
  T[, fragment := grepl("frag|partial|pseudo|relic", tip, ignore.case = TRUE) |
                  (!is.na(len) & is.finite(medlen) & len < 0.5 * medlen)]
  T[, `:=`(seqid = NA_character_, start = NA_integer_, end = NA_integer_)]
  for (s in intersect(unique(T$k), names(FULLNAME))) {
    g <- gff_genes(s)
    if (!nrow(g)) next
    for (i in which(T$k == s & !is.na(T$gene))) {
      hit <- g[id == T$gene[i]]
      if (!nrow(hit)) hit <- g[id_bare == sub("_(rnasupp|norna)$", "", T$gene[i])]
      if (nrow(hit)) set(T, i, c("seqid", "start", "end"), list(hit$seqid[1], hit$start[1], hit$end[1]))
    }
  }
  fwrite(T, file.path(OUT, sprintf("tips_%s.tsv", tag)), sep = "\t")

  hr(sprintf("tree: %s", tf))
  cat(sprintf("   %d tips | %s | alignment: %s | median Droseraceae length %s aa\n", length(labs),
              if (family) "alpha-kleisin FAMILY tree (holds SYN1-3 too)" else "SYN4-only tree",
              if (is.null(aln)) "none found beside the tree" else sprintf("%s (%.0f%% of tips)", basename(aln$file), 100 * aln$cov),
              if (is.finite(medlen)) format(medlen) else "?"))
  unk <- T[group == "unknown"]$tip
  if (length(unk)) cat(sprintf("   !! species not recognised for %d tips: %s\n", length(unk), paste(head(unk, 8), collapse = " ")))

  if (family) {
    # in the family tree SYN1-3 genes are expected in other OGs; show where the Drosera tips sit instead
    comp <- T[group %in% DROS, .N, by = .(og = ifelse(is.na(og), "no gene in run", og))][order(-N)]
    cat(sprintf("   Drosera tips by OG: %s\n", paste(sprintf("%s %d", comp$og, comp$N), collapse = ", ")))
  }

  # ------------------------------------------------ per Drosera species
  S <- rbindlist(lapply(intersect(DROS, unique(T$group)), function(s) {
    t <- T[group == s]
    if (family) t <- t[in_og == TRUE | is.na(og) | unassigned == TRUE]
    full <- t[fragment == FALSE]; frag <- t[fragment == TRUE]
    genes <- unique(na.omit(full$gene))
    g_in  <- unique(full[in_og == TRUE]$gene)
    oth   <- full[!is.na(og) & og != OGID & unassigned == FALSE]
    g_oth <- if (nrow(oth)) unique(paste0(oth$gene, "(", oth$og, ")")) else character(0)
    g_una <- unique(full[unassigned == TRUE]$gene)
    unmatched <- full[is.na(gene)]$tip
    missing <- setdiff(og_members[sp == s]$gene, unique(na.omit(t$gene)))
    loc <- unique(full[!is.na(seqid), .(gene, seqid)])
    if (!nrow(loc) && length(genes)) loc <- data.table(gene = character(), seqid = character())
    v <- if (s %in% names(violin)) violin[[s]] else NA_integer_
    iso <- full[!is.na(gene), .(n = .N, iso = paste(sort(isoform), collapse = "")), by = gene][n > 1]
    why <- character()
    if (nrow(iso)) why <- c(why, sprintf("%d extra tip(s) are isoforms of one gene [%s]", sum(iso$n - 1),
                                         paste(sprintf("%s:%s", iso$gene, iso$iso), collapse = " ")))
    if (length(g_oth)) why <- c(why, sprintf("%d gene(s) in other OGs [%s]", length(g_oth), paste(g_oth, collapse = " ")))
    if (length(g_una)) why <- c(why, sprintf("%d unassigned gene(s) [%s]", length(g_una), paste(g_una, collapse = " ")))
    if (length(unmatched)) why <- c(why, sprintf("%d tip(s) match no gene in this run [%s]", length(unmatched), paste(unmatched, collapse = " ")))
    if (length(missing)) why <- c(why, sprintf("%d %s gene(s) not in the tree [%s]", length(missing), OGID, paste(missing, collapse = " ")))
    if (nrow(frag)) why <- c(why, sprintf("%d fragment(s) [%s]", nrow(frag), paste(sprintf("%s=%s%s", frag$tip,
                                         ifelse(is.na(frag$og), "no gene in run", ifelse(frag$in_og, "IN", ifelse(frag$unassigned, "unassigned", frag$og))),
                                         ifelse(is.na(frag$len), "", sprintf(",%daa", frag$len))), collapse = " ")))
    agree <- !is.na(v) && length(g_in) == v && !length(g_oth) && !length(g_una) && !length(unmatched)
    data.table(species = s, cent = ifelse(s %in% HOLO, "holo", "mono"), violin = v,
               tips = nrow(t), full_tips = nrow(full), fragments = nrow(frag),
               genes = length(genes), in_OG = length(g_in), other_OG = length(g_oth), unassigned = length(g_una),
               unmatched = length(unmatched), OG_missing = length(missing),
               loci = loci_proxy(loc$seqid), chromosomes = paste(sort(unique(loc$seqid)), collapse = ","),
               verdict = paste0(if (agree) "AGREES" else "DIFFERS", if (length(why)) paste0(": ", paste(why, collapse = "; ")) else ""))
  }))
  if (nrow(S)) {
    fwrite(S, file.path(OUT, sprintf("summary_%s.tsv", tag)), sep = "\t")
    cat("   per species: tips = all; full = not a fragment; genes = distinct genes among full tips;\n")
    cat(sprintf("   in_OG = of those, genes in %s; loci = copies per haploid set (hap1/hap2 of one chromosome = 1)\n", OGID))
    print(S[, .(species, cent, violin, tips, full = full_tips, frag = fragments, genes, in_OG, other_OG,
                unassigned, unmatched, OG_missing, loci)], row.names = FALSE)
    cat("\n")
    for (i in seq_len(nrow(S))) cat(sprintf("   %-11s %s\n", S$species[i], S$verdict[i]))
    cat("   chromosomes carrying the full-length genes:\n")
    for (i in seq_len(nrow(S))) cat(sprintf("     %-11s %s\n", S$species[i], S$chromosomes[i]))
  }
  O <- T[!(group %in% DROS), .(tips = .N, genes = uniqueN(gene), frag = sum(fragment)), by = group]
  if (nrow(O)) cat(sprintf("   non-Drosera: %s\n", paste(sprintf("%s %d tips/%d genes%s", O$group, O$tips, O$genes,
                                                       ifelse(O$frag > 0, sprintf(" (%d fragment)", O$frag), "")), collapse = "; ")))

  # ------------------------------------------------ the daughter clades of the Droseraceae node
  outg <- which(T$group %in% OUTGROUPS)
  dros <- which(T$group %in% c(DROS, "Dionaea"))
  rooted <- NULL
  if (!family && length(outg) >= 1 && length(dros) >= 3) {
    ok <- tryCatch(length(outg) == 1 || is.monophyletic(tr, outg, reroot = TRUE), error = function(e) FALSE)
    if (ok) rooted <- tryCatch(root(tr, outgroup = outg, resolve.root = TRUE, edgelabel = TRUE), error = function(e) NULL)
    if (is.null(rooted)) {
      cat("   outgroup tips are not one clade: tree not rooted, no clade split\n")
    } else if (!is.monophyletic(rooted, tr$tip.label[dros])) {
      cat("   Droseraceae tips are not one clade once rooted on the outgroups: no clade split\n")
    } else {
      node <- getMRCA(rooted, tr$tip.label[dros])
      kids <- rooted$edge[rooted$edge[, 1] == node, 2]
      cat(sprintf("   Droseraceae node, rooted on %s, has %d daughter clades (distinct full-length genes; +frag = fragments):\n",
                  paste(unique(T$group[outg]), collapse = "+"), length(kids)))
      for (j in seq_along(kids)) {
        kt <- if (kids[j] <= Ntip(rooted)) rooted$tip.label[kids[j]] else extract.clade(rooted, kids[j])$tip.label
        K <- T[tip %in% kt]
        per <- K[, .(g = uniqueN(gene[!fragment & !is.na(gene)]) + sum(!fragment & is.na(gene)), f = sum(fragment)),
                 by = group][order(match(group, c(DROS, "Dionaea")))]
        sup <- if (kids[j] > Ntip(rooted) && !is.null(rooted$node.label)) rooted$node.label[kids[j] - Ntip(rooted)] else ""
        cat(sprintf("     clade %d (%2d tips%s): %s\n", j, length(kt), if (nzchar(sup)) paste0(", support ", sup) else "",
                    paste(sprintf("%s %d%s", per$group, per$g, ifelse(per$f > 0, sprintf("+%dfrag", per$f), "")), collapse = ", ")))
      }
    }
  }
  # compact topology to read by eye: species|tip id|OG status|FRAG
  rl <- tr
  rl$tip.label <- T[, paste0(group, "|", core, ifelse(in_og, "|IN", ifelse(is.na(og), "", paste0("|", og))), ifelse(fragment, "|FRAG", ""))]
  rl$node.label <- NULL
  if (!is.null(rooted)) rl <- tryCatch(root(rl, outgroup = outg, resolve.root = TRUE), error = function(e) rl)
  writeLines(write.tree(rl), file.path(OUT, sprintf("topology_%s.nwk", tag)))
}
cat(sprintf("\ntables: %s/{tips,summary,topology}_*\n", OUT))
