#!/usr/bin/env Rscript
# syn4_tree_vs_og.R (v2) -- cut the SYN4 clade out of the alpha-kleisin family tree and count it three ways
#
# Why v2: the orthogroups on disk come from the 12 Aug 2026 rerun (11 species, tagged annotation). OrthoFinder
# renumbers orthogroups on every run, so "OG0001569" in that run is an unrelated family; the run the violin
# used (8 Drosera, before the rerun) no longer exists. So the clade is compared against:
#   fig      the violin's numbers, typed in from the figure (FIG_COUNTS)
#   old      the old OG0001569 members that survive in a file written before the rerun (OLD_MEMBERS)
#   current  the orthogroup of today's run that holds the clade's Drosera genes (this one includes scorpioides)
# Per tip it reports species, gene, isoform, current orthogroup, chromosome, length and fragment status, so a
# mismatch reads as isoforms, genes split into other orthogroups, genes missing from the tree, or fragments.
#
# SYN4 clade = walk up from the SYN4 query tip(s) (AT_SYN4 / AT5G16270) in the tree rooted on SYN1, and keep
# the largest clade that still holds no SYN1/SYN2/SYN3 marker tip.
#
# Usage: Rscript syn4_tree_vs_og.R OUT_DIR ORTHOGROUPS_DIR REPO TREE FIG_COUNTS [OLD_MEMBERS]
#   FIG_COUNTS  e.g. regia=1,paradoxa=2,roseana=2,filiformis=3,binata=4,capensis=4,aliciae=6,tokaiensis=6
suppressPackageStartupMessages({ library(ape); library(data.table) })
options(width = 250)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("usage: OUT_DIR ORTHOGROUPS_DIR REPO TREE FIG_COUNTS [OLD_MEMBERS]")
OUT <- args[1]; OGD <- args[2]; REPO <- args[3]; TF <- args[4]; FIGS <- args[5]
OLDF <- if (length(args) >= 6) args[6] else NA_character_
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

DROS  <- c("regia", "paradoxa", "roseana", "scorpioides", "binata", "capensis", "aliciae", "tokaiensis", "filiformis")
HOLO  <- c("regia", "paradoxa", "roseana", "scorpioides")
OTHER <- c(arabidopsis = "Arabidopsis", atha = "Arabidopsis", ath = "Arabidopsis", athal = "Arabidopsis", at = "Arabidopsis",
           beta = "Beta", bvul = "Beta", bvulgaris = "Beta",
           nepenthes = "Nepenthes", ngra = "Nepenthes", ngracilis = "Nepenthes",
           dionaea = "Dionaea", dmus = "Dionaea", dmuscipula = "Dionaea")
FULLNAME <- c(setNames(paste0("Drosera_", DROS), DROS), dionaea = "Dionaea_muscipula", nepenthes = "Nepenthes_gracilis")
SYN4_SEED <- "(^|[_|])SYN4$|AT5G16270"
OTHER_KLEISIN <- "(^|[_|])SYN[123]$|AT5G05490|AT5G40840|AT3G59550"
SYN1_ROOT <- "(^|[_|])SYN1$|AT5G05490"
hr <- function(t) cat(sprintf("\n== %s\n", t))
short_sp <- function(col) { s <- tolower(sub("^Drosera_", "", col)); s[s == "dionaea_muscipula"] <- "dionaea"; s[s == "nepenthes_gracilis"] <- "nepenthes"; s }

fig <- local({
  kv <- strsplit(strsplit(FIGS, ",", fixed = TRUE)[[1]], "=", fixed = TRUE)
  setNames(as.integer(vapply(kv, `[`, "", 2)), vapply(kv, `[`, "", 1))
})

# ---------------------------------------------------------------- today's orthogroups: gene -> OG
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
og[, sp := short_sp(col)]
og[, gene_bare := sub("_(rnasupp|norna)$", "", gene)]
by_gene <- og[, .(key = paste(sp, gene), gene, Orthogroup, unassigned)]
by_bare <- og[, .(key = paste(sp, gene_bare), gene, Orthogroup, unassigned)]
setkey(by_gene, key); setkey(by_bare, key)
run_species <- sort(unique(og$sp))
gcount <- fread(file.path(OGD, "Orthogroups.GeneCount.tsv"), sep = "\t", header = TRUE)
setnames(gcount, names(gcount), c("Orthogroup", short_sp(names(gcount)[-1])))
hr("today's OrthoFinder run (the one on disk)")
cat(sprintf("   %s\n   species (%d): %s\n   %s orthogroups; Orthogroups.tsv modified %s\n", OGD, length(run_species),
            paste(run_species, collapse = " "), format(uniqueN(og[unassigned == FALSE]$Orthogroup), big.mark = ","),
            format(file.mtime(ogf), "%F %H:%M")))

# ---------------------------------------------------------------- helpers
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
  if (is.na(k) || !(k %in% run_species) || is.na(core) || core == "") return(none)
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
  # copies per haploid set: chrN_hap1 + chrN_hap2 = one locus seen twice; h1tg/h2tg contigs can't be paired
  # by name, so take the larger haplotype's count; collapsed assemblies: every copy counts
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
clade_tips <- function(tr, node) if (node <= Ntip(tr)) tr$tip.label[node] else extract.clade(tr, node)$tip.label

# ---------------------------------------------------------------- old OG0001569 members (written before the rerun)
old <- NULL
if (!is.na(OLDF) && file.exists(OLDF)) {
  d <- fread(OLDF, colClasses = "character")
  pick_sp <- function(v) { v <- tolower(v); out <- rep(NA_character_, length(v))
    for (k in DROS) out[is.na(out) & grepl(k, v, fixed = TRUE)] <- k; out }
  spcol <- names(d)[which.max(vapply(d, function(v) mean(!is.na(pick_sp(v))), 0))]
  d[, sp_ := pick_sp(get(spcol))]
  best <- NA_character_; best_n <- 0L
  for (cn in setdiff(names(d), c(spcol, "sp_"))) {
    n <- sum(vapply(seq_len(nrow(d)), function(i) !is.na(match_gene(d$sp_[i], strip_prefix(d[[cn]][i], d$sp_[i]))$gene), TRUE))
    if (n > best_n) { best_n <- n; best <- cn }
  }
  hr(sprintf("old OG0001569 members: %s", OLDF))
  cat(sprintf("   %d rows; species column '%s'; gene column %s\n", nrow(d), spcol,
              if (is.na(best)) "none found (no column matches today's gene IDs)" else sprintf("'%s' (%d of %d rows matched to today's genes)", best, best_n, nrow(d))))
  old <- d[, .(sp = sp_, old_id = if (is.na(best)) NA_character_ else get(best))]
  old[, gene := vapply(seq_len(.N), function(i) if (is.na(old_id[i])) NA_character_ else match_gene(sp[i], strip_prefix(old_id[i], sp[i]))$gene, "")]
  old_keys <- paste(old$sp, old$gene)                       # built outside by_gene[...]: by_gene has its own 'gene' column
  old$og_today <- by_gene[J(old_keys), Orthogroup, mult = "first"]
  old[is.na(gene), og_today := NA_character_]
  cmp <- merge(data.table(sp = names(fig), fig = as.integer(fig)), old[, .(old_rows = .N), by = sp], by = "sp", all = TRUE)
  cat("   rows per species vs the violin figure:", paste(sprintf("%s %s/%s", cmp$sp, cmp$old_rows, cmp$fig), collapse = ", "), "\n")
  if (!is.na(best)) cat("   where those genes sit today:", paste(old[, .N, by = .(og = ifelse(is.na(og_today), "not found", og_today))][order(-N), sprintf("%s %d", og, N)], collapse = ", "), "\n")
}

# ---------------------------------------------------------------- the tree
tr <- read.tree(TF)
if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
labs <- tr$tip.label
aln <- find_alignment(TF, labs)
T <- data.table(tip = labs)
T[, k := vapply(tip, species_of, "")]
T[, group := ifelse(k %in% DROS, k, ifelse(is.na(k), "unknown", OTHER[k]))]
T[, core := unname(mapply(strip_prefix, tip, k))]
T[, isoform := ifelse(grepl("\\.t[0-9]+", core), sub("^.*(\\.t[0-9]+).*$", "\\1", core), "")]
m <- lapply(seq_len(nrow(T)), function(i) match_gene(T$k[i], T$core[i]))
T[, gene := vapply(m, function(z) z$gene, "")]
T[, og := vapply(m, function(z) z$og, "")]
T[, unassigned := vapply(m, function(z) isTRUE(z$unassigned), TRUE)]
T[is.na(gene) & group %in% c("Arabidopsis", "Beta", "Nepenthes", "Dionaea"), gene := sub("\\.[0-9]+$", "", sub("\\.t[0-9]+$", "", core))]
T[, len := NA_integer_]
if (!is.null(aln)) T[, len := nchar(aln$fa$seq[match(norm(tip), aln$fa$nid)])]
medlen <- suppressWarnings(median(T[group %in% c(DROS, "Dionaea") & !is.na(len)]$len))
T[, fragment := grepl("frag|partial|pseudo|relic", tip, ignore.case = TRUE) | (!is.na(len) & is.finite(medlen) & len < 0.5 * medlen)]
T[, `:=`(seqid = NA_character_, start = NA_integer_)]
for (s in intersect(unique(T$k), names(FULLNAME))) {
  g <- gff_genes(s)
  if (!nrow(g)) next
  for (i in which(T$k == s & !is.na(T$gene))) {
    hit <- g[id == T$gene[i]]
    if (!nrow(hit)) hit <- g[id_bare == sub("_(rnasupp|norna)$", "", T$gene[i])]
    if (nrow(hit)) set(T, i, c("seqid", "start"), list(hit$seqid[1], hit$start[1]))
  }
}
hr(sprintf("tree: %s", TF))
cat(sprintf("   %d tips; alignment %s; median Droseraceae length %s aa\n", length(labs),
            if (is.null(aln)) "not found" else sprintf("%s (%.0f%% of tips)", basename(aln$file), 100 * aln$cov), format(medlen)))
unk <- T[group == "unknown"]$tip
if (length(unk)) cat(sprintf("   !! species not recognised for %d tips: %s\n", length(unk), paste(head(unk, 8), collapse = " ")))

# ---------------------------------------------------------------- root on SYN1, cut out the SYN4 clade
r1 <- grep(SYN1_ROOT, labs, ignore.case = TRUE)
seed <- grep(SYN4_SEED, labs, ignore.case = TRUE)
if (!length(seed)) stop("no SYN4 query tip (AT_SYN4 / AT5G16270) in this tree")
rt <- if (length(r1)) root(tr, outgroup = r1[1], resolve.root = TRUE, edgelabel = TRUE) else tr
cat(sprintf("   rooted on %s; SYN4 seed tip(s): %s\n", if (length(r1)) labs[r1[1]] else "nothing (no SYN1 tip)",
            paste(labs[seed], collapse = " ")))
seed_lab <- labs[seed]
node <- if (length(seed_lab) > 1) getMRCA(rt, seed_lab) else rt$edge[rt$edge[, 2] == which(rt$tip.label == seed_lab), 1]
bad <- function(nd) any(grepl(OTHER_KLEISIN, clade_tips(rt, nd), ignore.case = TRUE))
if (bad(node)) stop("the smallest clade around the SYN4 seed already holds SYN1/2/3 markers -- check the tree")
repeat {
  up <- rt$edge[rt$edge[, 2] == node, 1]
  if (!length(up) || bad(up)) break
  node <- up
}
S4 <- clade_tips(rt, node)
sup <- function(nd) if (nd > Ntip(rt) && !is.null(rt$node.label)) rt$node.label[nd - Ntip(rt)] else ""
T[, in_syn4 := tip %in% S4]
fwrite(T, file.path(OUT, "tips_syn4_clade.tsv"), sep = "\t")
describe <- function(tips) {
  K <- T[tip %in% tips]
  per <- K[, .(g = uniqueN(gene[!fragment & !is.na(gene)]) + sum(!fragment & is.na(gene)), f = sum(fragment)), by = group]
  per <- per[order(match(group, c(DROS, "Dionaea", "Nepenthes", "Beta", "Arabidopsis", "unknown")))]
  paste(sprintf("%s %d%s", per$group, per$g, ifelse(per$f > 0, sprintf("+%dfrag", per$f), "")), collapse = ", ")
}
hr(sprintf("SYN4 clade: %d tips (support %s)", length(S4), sup(node)))
kids <- rt$edge[rt$edge[, 1] == node, 2]
for (j in seq_along(kids)) cat(sprintf("   daughter %d (%d tips, support %s): %s\n", j, length(clade_tips(rt, kids[j])), sup(kids[j]), describe(clade_tips(rt, kids[j]))))
dr <- T[in_syn4 & group %in% c(DROS, "Dionaea")]$tip
if (length(dr) >= 2 && is.monophyletic(rt, dr)) {
  dn <- getMRCA(rt, dr)
  cat(sprintf("   Droseraceae part is one clade (support %s); its daughters:\n", sup(dn)))
  for (kd in rt$edge[rt$edge[, 1] == dn, 2]) cat(sprintf("     %d tips, support %s: %s\n", length(clade_tips(rt, kd)), sup(kd), describe(clade_tips(rt, kd))))
} else cat("   Droseraceae tips in the SYN4 clade are NOT one clade (outgroup sequences sit among them)\n")
cat("   (numbers = distinct full-length genes; +frag = fragments)\n")

# ---------------------------------------------------------------- the clade's Drosera tips, one line each
D <- T[in_syn4 == TRUE & group %in% DROS]
cur_og <- names(sort(table(D[fragment == FALSE & !is.na(og) & !unassigned]$og), decreasing = TRUE))[1]
hr(sprintf("Drosera tips in the SYN4 clade (today's orthogroup holding most of them: %s)", cur_og))
setorder(D, group, gene, isoform)
for (i in seq_len(nrow(D))) cat(sprintf("   %-10s %-38s %-26s %-12s %-16s %5s aa%s\n", D$group[i], D$tip[i],
  ifelse(is.na(D$gene[i]), "-", D$gene[i]), ifelse(is.na(D$og[i]), "no gene", ifelse(D$unassigned[i], "unassigned", D$og[i])),
  ifelse(is.na(D$seqid[i]), "-", D$seqid[i]), ifelse(is.na(D$len[i]), "?", D$len[i]), ifelse(D$fragment[i], "  FRAGMENT", "")))

# ---------------------------------------------------------------- per species: fig vs old vs tree vs today
S <- rbindlist(lapply(DROS, function(s) {
  t <- D[group == s]; full <- t[fragment == FALSE]
  genes <- unique(na.omit(full$gene))
  o <- if (!is.null(old)) old[sp == s] else data.table(old_id = character(), gene = character())
  og_today <- if (!is.na(cur_og) && s %in% names(gcount)) as.integer(gcount[Orthogroup == cur_og][[s]]) else NA_integer_
  cur_genes <- og[Orthogroup == cur_og & sp == s & unassigned == FALSE]$gene
  iso <- full[!is.na(gene), .(n = .N, iso = paste(sort(isoform), collapse = "")), by = gene][n > 1]
  oth <- full[!is.na(og) & og != cur_og]
  why <- character()
  if (nrow(iso)) why <- c(why, sprintf("isoform tips: %s", paste(sprintf("%s:%s", iso$gene, iso$iso), collapse = " ")))
  if (nrow(oth)) why <- c(why, sprintf("in other orthogroups today: %s", paste(unique(sprintf("%s(%s)", oth$gene, ifelse(oth$unassigned, "unassigned", oth$og))), collapse = " ")))
  if (nrow(o) && any(!is.na(o$gene))) {
    om <- setdiff(na.omit(o$gene), genes); tn <- setdiff(genes, na.omit(o$gene))
    if (length(om)) why <- c(why, sprintf("old members not in the clade: %s", paste(om, collapse = " ")))
    if (length(tn)) why <- c(why, sprintf("clade genes not among old members: %s", paste(tn, collapse = " ")))
  }
  cm <- setdiff(cur_genes, unique(na.omit(t$gene)))
  if (length(cm)) why <- c(why, sprintf("%s genes not in the tree: %s", cur_og, paste(cm, collapse = " ")))
  if (any(t$fragment)) why <- c(why, sprintf("fragments: %s", paste(sprintf("%s(%saa)", t[fragment == TRUE]$tip, t[fragment == TRUE]$len), collapse = " ")))
  data.table(species = s, cent = ifelse(s %in% HOLO, "holo", "mono"),
             fig = if (s %in% names(fig)) fig[[s]] else NA_integer_,
             old = if (nrow(o)) nrow(o) else NA_integer_,
             tips = nrow(t), frag = sum(t$fragment), genes = length(genes),
             loci = loci_proxy(unique(full[!is.na(seqid), .(gene, seqid)])$seqid),
             today_OG = og_today, notes = paste(why, collapse = "; "))
}))
fwrite(S, file.path(OUT, "summary_syn4_clade.tsv"), sep = "\t")
hr("per species. fig = violin figure; old = old OG0001569 rows; tips/frag/genes/loci = SYN4 clade of this tree;")
cat(sprintf("   today_OG = genes in %s, today's orthogroup for this clade (scorpioides has no tips in the tree)\n", cur_og))
print(S[, .(species, cent, fig, old, tips, frag, genes, loci, today_OG)], row.names = FALSE)
cat("\n")
for (i in which(nzchar(S$notes))) cat(sprintf("   %-11s %s\n", S$species[i], S$notes[i]))
if (!is.na(cur_og)) {
  r <- gcount[Orthogroup == cur_og]
  cat(sprintf("\n   %s in today's run, all species: %s\n", cur_og,
              paste(sprintf("%s %s", names(r)[-1], unlist(r[, -1])), collapse = ", ")))
}
rl <- keep.tip(rt, S4)
rl$tip.label <- T[match(rl$tip.label, tip), paste0(group, "|", core, ifelse(is.na(og), "", paste0("|", og)), ifelse(fragment, "|FRAG", ""))]
writeLines(write.tree(rl), file.path(OUT, "topology_syn4_clade.nwk"))
cat(sprintf("\ntables: %s/{tips,summary,topology}_syn4_clade.*\n", OUT))
