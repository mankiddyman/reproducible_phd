#!/usr/bin/env Rscript
# tpxl6_first_exon.R -- the two R steps of tpxl6_first_exon.sh (what the whole thing does is described there)
#   locate REPO SPECIES GENOME OUTDIR BINDIR THREADS
#            for every TPXL6 copy aligned by miniprot from after amino acid 1: align the start of the complete copy
#            (3 kb upstream of its first exon to the end of exon 2, CARP repeats masked) to the 200 kb upstream of the
#            incomplete copy (minimap2 -x asm20), and scan the same DNA base by base, both strands, for the complete
#            copy's first exon (E = expected number of matches this good by chance in that DNA, from its base
#            composition). Writes OUTDIR/anchor.tsv (one row per incomplete copy), blocks.tsv, exon1_hits.tsv, and
#            rescan_region.tsv + usable_bams.tsv for the shell step
#   report REPO SPECIES GENOME OUTDIR BINDIR THREADS OUT_PREFIX
#            CpG methylation around the found start, the complete copy's start and the related gene, in the windows of
#            tpxl6_methyl_tss_check.R (CpGs inside CARP TEs left out); OUT_PREFIX.tsv/.png/.pdf
# Coordinates: "x" is bp along a gene, 5'->3', from its first coding base as miniprot found it (for the incomplete
# copy that is exon 2); the found start is given as x on the incomplete copy (negative = upstream of its exon 2).
# FIG_BANNER (environment variable): when set, printed in red at the top right of the figure, e.g. to mark test data.
suppressPackageStartupMessages(library(data.table)); options(width = 220, scipen = 20)
a <- commandArgs(TRUE)
if (length(a) < 7 || !a[1] %in% c("locate", "report"))
  stop("usage: tpxl6_first_exon.R locate|report REPO SPECIES GENOME OUTDIR BINDIR THREADS [OUT_PREFIX]")
MODE <- a[1]; REPO <- a[2]; SP <- a[3]; G <- a[4]; O <- a[5]; BIN <- a[6]; TH <- a[7]; OUTP <- if (length(a) >= 8) a[8] else ""
if (!startsWith(G, "/")) G <- file.path(REPO, G)
SAM <- file.path(BIN, "samtools"); MM2 <- file.path(BIN, "minimap2")
H  <- file.path(REPO, "results/comparative/holocentricity")
MP <- file.path(H, "tpx2/miniprot"); MD <- file.path(H, "tpx2/methyl"); D <- file.path(MD, SP)
CARP <- file.path(REPO, "results", SP, "repeats/carp/Repeat_Annotation_Unified.gff3")
Q6 <- "regia_g4881.t1"; Q5 <- "regia_g4900.t1"
UP <- 3000; SEARCH <- 200000; FL <- 5000; MAXMM <- 10
X0 <- -3000; X1 <- 1500; SMOOTH <- 100
WIN <- list(`-2000..-1000` = c(-2000, -1000), `-1000..-500` = c(-1000, -500), `-500..0` = c(-500, 0), `0..+500` = c(0, 500))
INK <- "#1f1f1e"; INK2 <- "#52514e"; MUTED <- "#8a8984"; RED <- "#b2182b"; RIB <- "#8a8984"
GCOL <- c(LTR = "#2a78d6", oTE = "#e87ba4", DNA = "#eda100", rep = "#cfcec9")

alpha <- function(col, a) { m <- col2rgb(col) / 255; rgb(m[1], m[2], m[3], a) }
nice <- function(sp) sub("^Drosera_", "D. ", sp)
hap  <- function(ctg) { h <- regmatches(ctg, regexpr("hap[0-9]+|h[0-9]+tg", ctg)); if (length(h)) sub("^h([0-9]+)tg$", "hap\\1", h) else ctg }
attr1 <- function(x, k) ifelse(grepl(paste0("(^|;)", k, "="), x), sub(paste0(".*(^|;)", k, "=([^;]*).*"), "\\2", x), NA)
first_of <- function(x, y) { x <- as.character(x); y <- as.character(y); x[is.na(x)] <- y[is.na(x)]; x }
kb <- function(x) sprintf("%.1f kb", x / 1000)
GFFCOLS <- c("contig", "src", "type", "start", "end", "score", "strand", "phase", "attr")
gff_text <- function(x) { x <- x[!startsWith(x, "#") & nzchar(x)]; if (!length(x)) return(NULL)
  fread(text = paste0(paste(x, collapse = "\n"), "\n"), header = FALSE, sep = "\t", quote = "", col.names = GFFCOLS) }   # one line alone would be read as a file name
grp_of <- function(cl) fifelse(startsWith(cl, "Class_I/LTR"), "LTR", fifelse(startsWith(cl, "Class_II"), "DNA",
                      fifelse(startsWith(cl, "Class_I"), "oTE",
                      fifelse(startsWith(cl, "Simple_repeat") | startsWith(cl, "Low_complexity"), NA_character_, "rep"))))
inside <- function(p, s, e) { o <- rep(FALSE, length(p)); for (i in seq_along(s)) o <- o | (p >= s[i] & p < e[i]); o }
rc <- function(s) paste(rev(strsplit(chartr("ACGTN", "TGCAN", s), "")[[1]]), collapse = "")
B4 <- c("T", "C", "A", "G")
CODON <- setNames(strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]],
                  paste0(rep(B4, each = 16), rep(rep(B4, each = 4), 4), rep(B4, 16)))
translate <- function(s) { if (is.na(s)) return(NA_character_); n <- nchar(s) %/% 3; if (!n) return("")
  x <- CODON[substring(s, seq(1, by = 3, length.out = n), seq(3, by = 3, length.out = n))]; x[is.na(x)] <- "X"; paste(x, collapse = "") }

# genome access; a "frame" is list(contig, strand, g0): x = 0 at genomic g0, x grows 5'->3' along the gene
fai <- fread(paste0(G, ".fai"), header = FALSE, select = 1:2, col.names = c("contig", "len"))
clen <- function(ctg) fai$len[fai$contig == ctg]
faidx <- function(ctg, lo, hi) { x <- system2(SAM, c("faidx", shQuote(G), shQuote(sprintf("%s:%d-%d", ctg, as.integer(lo), as.integer(hi)))), stdout = TRUE)
  toupper(paste(x[!startsWith(x, ">")], collapse = "")) }
tox <- function(fr, p) if (fr$strand == "+") p - fr$g0 else fr$g0 - p
tog <- function(fr, x) if (fr$strand == "+") fr$g0 + x else fr$g0 - x
grange <- function(fr, from, to) sort(c(tog(fr, from), tog(fr, to - 1)))          # genomic [lo, hi] of x in [from, to)
gseq <- function(fr, from, to) {                                                    # sequence of x in [from, to), clipped to the contig
  L <- clen(fr$contig)
  if (fr$strand == "+") { lo <- max(1, fr$g0 + from); hi <- min(L, fr$g0 + to - 1)
    list(seq = faidx(fr$contig, lo, hi), from = lo - fr$g0, to = hi - fr$g0 + 1)
  } else { lo <- max(1, fr$g0 - (to - 1)); hi <- min(L, fr$g0 - from)
    list(seq = rc(faidx(fr$contig, lo, hi)), from = fr$g0 - hi, to = fr$g0 - lo + 1) }
}
# CARP Level-1 features (no Parent=) overlapping any of the genomic ranges, one pass over the file
carp_get <- function(rng) {
  if (!file.exists(CARP) || !nrow(rng)) return(NULL)
  rf <- tempfile(fileext = ".tsv"); fwrite(rng[, .(contig, as.integer(lo), as.integer(hi))], rf, sep = "\t", col.names = FALSE)
  x <- system(sprintf("awk -F'\\t' -v rf=%s 'FILENAME == rf { n++; C[n] = $1; L[n] = $2; H[n] = $3; next }
    /^#/ || $9 ~ /(^|;)Parent=/ { next } { for (i = 1; i <= n; i++) if ($1 == C[i] && $4 <= H[i] && $5 >= L[i]) { print; break } }' %s %s",
    shQuote(rf), shQuote(rf), shQuote(CARP)), intern = TRUE)
  unlink(rf)
  d <- gff_text(x); if (is.null(d)) return(NULL)
  d[, cl := first_of(attr1(attr, "classification"), attr1(attr, "Name"))][, grp := grp_of(cl)]
  d
}
bed_read <- function(f) {
  if (is.na(f) || !file.exists(f)) return(NULL)
  b <- fread(cmd = paste("gzip -dc", shQuote(f), "| grep -v '^#'"), header = FALSE, sep = "\t")
  if (!nrow(b)) return(NULL)
  setnames(b, 1:6, c("contig", "begin", "end", "mod", "type", "cov")); b[, pos := begin + 1L][, .(contig, pos, mod, cov)]
}
paf_read <- function(x) {
  if (!length(x)) return(data.table())
  rbindlist(lapply(strsplit(x, "\t"), function(f) data.table(qs = as.numeric(f[3]), qe = as.numeric(f[4]), str = f[5], tname = f[6],
    ts = as.numeric(f[8]), te = as.numeric(f[9]), nm = as.numeric(f[10]), bl = as.numeric(f[11]), mapq = as.integer(f[12]),
    tp = sub("^tp:A:", "", grep("^tp:A:", f, value = TRUE)[1]), cs = sub("^cs:Z:", "", grep("^cs:Z:", f, value = TRUE)[1]))))
}
# target offset (0-based) aligned to query offset qpos (0-based) in a + strand alignment, by walking its cs string
cs_map <- function(cs, qs, ts, qpos) {
  ops <- regmatches(cs, gregexpr(":[0-9]+|\\*[a-z][a-z]|\\+[a-z]+|-[a-z]+", cs))[[1]]
  q <- qs; t <- ts
  for (op in ops) {
    k <- substr(op, 1, 1)
    if (k == ":") { n <- as.integer(substring(op, 2)); if (qpos < q + n) return(t + qpos - q); q <- q + n; t <- t + n }
    else if (k == "*") { if (qpos == q) return(t); q <- q + 1; t <- t + 1 }
    else if (k == "+") { n <- nchar(op) - 1; if (qpos < q + n) return(t); q <- q + n }
    else if (k == "-") { t <- t + nchar(op) - 1 }
  }
  NA_real_
}
# a + strand alignment cut at every gap >= MINGAP bp into segments (0-based query/target starts and ends, matches,
# mismatches), so an insertion such as a TE shows up as a gap between segments, not as one stretched block
MINGAP <- 100
cs_segments <- function(cs, qs, ts) {
  ops <- regmatches(cs, gregexpr(":[0-9]+|\\*[a-z][a-z]|\\+[a-z]+|-[a-z]+", cs))[[1]]
  q <- qs; t <- ts; sq <- qs; st <- ts; m <- 0; x <- 0; out <- list()
  for (op in ops) {
    k <- substr(op, 1, 1)
    if (k == ":") { n <- as.integer(substring(op, 2)); q <- q + n; t <- t + n; m <- m + n }
    else if (k == "*") { q <- q + 1; t <- t + 1; x <- x + 1 }
    else { n <- nchar(op) - 1
      if (n >= MINGAP) { if (m + x > 0) out[[length(out) + 1]] <- c(sq, q, st, t, m, x)
                         if (k == "+") q <- q + n else t <- t + n; sq <- q; st <- t; m <- 0; x <- 0 }
      else if (k == "+") q <- q + n else t <- t + n }
  }
  if (m + x > 0) out[[length(out) + 1]] <- c(sq, q, st, t, m, x)
  if (!length(out)) return(data.table(sqs = numeric(), sqe = numeric(), sts = numeric(), ste = numeric(), match = numeric(), mism = numeric()))
  setnames(as.data.table(do.call(rbind, out)), c("sqs", "sqe", "sts", "ste", "match", "mism"))
}
# % identity over aligned bases (gaps left out), from a cs string
cs_ident <- function(cs) { ops <- regmatches(cs, gregexpr(":[0-9]+|\\*[a-z][a-z]", cs))[[1]]; k <- substr(ops, 1, 1)
  m <- sum(as.integer(substring(ops[k == ":"], 2))); x <- sum(k == "*"); if (m + x > 0) round(100 * m / (m + x), 1) else NA_real_ }
# mismatches of probe p at every offset of sequence s (N counts as a mismatch)
ham <- function(s, p) {
  sv <- match(strsplit(s, "")[[1]], c("A", "C", "G", "T")); pv <- match(strsplit(p, "")[[1]], c("A", "C", "G", "T"))
  n <- length(sv) - length(pv) + 1; if (n < 1) return(integer())
  mm <- integer(n); for (i in seq_along(pv)) { v <- sv[i:(i + n - 1)]; mm <- mm + (is.na(v) | v != pv[i]) }
  mm
}

# TPXL6 copies (miniprot) and the related gene
L <- fread(file.path(MP, "tpxl6_miniprot_loci.tsv"))
g <- gff_text(readLines(file.path(MP, paste0(SP, ".gff"))))
copy_of <- function(r, q) {
  m <- g[type == "mRNA" & contig == r$contig & start == r$start & end == r$end & startsWith(attr1(attr, "Target"), paste0(q, " "))][1]
  cds <- g[type == "CDS" & attr1(attr, "Parent") == attr1(m$attr, "ID")]
  fr <- list(contig = r$contig, strand = r$strand, g0 = if (r$strand == "+") r$start else r$end)
  tg <- tstrsplit(attr1(cds$attr, "Target"), " ")
  ex <- data.table(x0 = pmin(tox(fr, cds$start), tox(fr, cds$end)), x1 = pmax(tox(fr, cds$start), tox(fr, cds$end)) + 1,
                   qs = as.integer(tg[[2]]), qe = as.integer(tg[[3]]))[order(x0)]
  c(fr, list(start = r$start, end = r$end, ex = ex, first_aa = min(ex$qs), score = r$score, lab = hap(r$contig)))
}
cps <- L[species == SP & query == Q6 & assigned == "TPXL6"]
if (!nrow(cps)) stop("no TPXL6 copy of ", SP, " in tpxl6_miniprot_loci.tsv")
CC <- lapply(seq_len(nrow(cps)), function(i) copy_of(cps[i], Q6))
FA <- sapply(CC, `[[`, "first_aa"); INC <- which(FA > 1); FULL <- which(FA == 1)
REF <- if (length(FULL)) CC[[FULL[which.max(sapply(CC[FULL], `[[`, "score"))]]] else NULL
chk <- fread(file.path(MD, "methyl_bam_check.tsv")); USE <- chk[species == SP & usable == TRUE]$bam_n
bed_name <- function(dir) if (length(USE) == 2) file.path(dir, paste0(SP, ".combined.bed.gz")) else
                          if (length(USE) == 1) file.path(dir, sprintf("%s.hap%d.bed.gz", SP, USE)) else NA_character_
seg_of <- function(ref) { e <- ref$ex; s <- data.table(part = "upstream", a = -Inf, b = e$x0[1])
  s <- rbind(s, data.table(part = "first exon", a = e$x0[1], b = e$x1[1]))
  if (nrow(e) > 1) s <- rbind(s, data.table(part = c("intron 1", "exon 2", "after exon 2"), a = c(e$x1[1], e$x0[2], e$x1[2]), b = c(e$x0[2], e$x1[2], Inf)))
  s }
parts_of <- function(segs, x0, x1) { o <- pmax(0, pmin(x1, segs$b) - pmax(x0, segs$a))
  paste(sprintf("%s %s bp", segs$part[o > 0], format(o[o > 0], big.mark = ",")), collapse = " + ") }

# ---------------------------------------------------------------------------------------------------------- locate
if (MODE == "locate") {
  cat(sprintf("   TPXL6 copies: %s\n", paste(sprintf("%s %s:%d-%d (%s), aligned from amino acid %d", sapply(CC, `[[`, "lab"),
      sapply(CC, `[[`, "contig"), sapply(CC, `[[`, "start"), sapply(CC, `[[`, "end"), sapply(CC, `[[`, "strand"), FA), collapse = "; ")))
  AN <- list(); BL <- list(); HI <- list(); RS <- list()
  if (!length(INC)) cat("   every copy is aligned from amino acid 1: nothing to find\n")
  if (length(INC) && is.null(REF)) stop("no copy of ", SP, " is aligned from amino acid 1, so there is no complete start to search with")
  bed <- bed_read(bed_name(D))
  for (k in INC) {
    tc <- CC[[k]]; ref <- REF; segs <- seg_of(ref)
    ex1 <- ref$ex[1]; L1 <- ex1$x1 - ex1$x0
    cat(sprintf("\n== %s (%s:%d-%d, %s strand) is aligned from amino acid %d; complete copy: %s (%s:%d-%d, %s strand), first exon %d bp (amino acids %d-%d)\n",
        tc$lab, tc$contig, tc$start, tc$end, tc$strand, tc$first_aa, ref$lab, ref$contig, ref$start, ref$end, ref$strand, L1, ex1$qs, ex1$qe))
    qto <- if (nrow(ref$ex) > 1) ref$ex$x1[2] + 200 else ex1$x1 + 2000
    q <- gseq(ref, ex1$x0 - UP, qto); t <- gseq(tc, -SEARCH, 2000)
    rq <- grange(ref, q$from, q$to); rt <- grange(tc, t$from, t$to)
    rp <- carp_get(data.table(contig = c(ref$contig, tc$contig), lo = c(rq[1], rt[1]), hi = c(rq[2], rt[2])))
    qch <- strsplit(q$seq, "")[[1]]                     # the complete copy's start, every CARP feature masked
    if (!is.null(rp)) { r2 <- rp[contig == ref$contig]
      for (j in seq_len(nrow(r2))) { x0 <- min(tox(ref, r2$start[j]), tox(ref, r2$end[j])); x1 <- max(tox(ref, r2$start[j]), tox(ref, r2$end[j])) + 1
        a0 <- max(1, x0 - q$from + 1); a1 <- min(length(qch), x1 - q$from); if (a0 <= a1) qch[a0:a1] <- "N" } }
    qf <- file.path(O, paste0("start_", ref$lab, ".fa")); tf <- file.path(O, paste0("upstream_", tc$lab, ".fa"))
    writeLines(c(paste0(">", ref$lab, "_start"), paste(qch, collapse = "")), qf)
    writeLines(c(paste0(">", tc$lab, "_upstream"), t$seq), tf)
    paf <- paf_read(system2(MM2, c("-cx", "asm20", "--cs", "-N", "50", "-p", "0.1", "-s", "100", "-t", TH, shQuote(tf), shQuote(qf)),
                            stdout = TRUE, stderr = FALSE))
    cat(sprintf("   masked %s of %s bp of the %s query (CARP repeats); searched %s:%d-%d (x = %d to %d on %s)\n",
        format(sum(qch == "N"), big.mark = ","), format(length(qch), big.mark = ","), ref$lab, tc$contig, rt[1], rt[2], t$from, t$to - 1, tc$lab))
    if (nrow(paf)) {
      sg <- rbindlist(lapply(seq_len(nrow(paf)), function(j) { p <- paf[j]
        s <- if (p$str == "+") cs_segments(p$cs, p$qs, p$ts) else {
          ops <- regmatches(p$cs, gregexpr(":[0-9]+|\\*[a-z][a-z]", p$cs))[[1]]; k <- substr(ops, 1, 1)
          data.table(sqs = p$qs, sqe = p$qe, sts = p$ts, ste = p$te, match = sum(as.integer(substring(ops[k == ":"], 2))), mism = sum(k == "*")) }
        if (nrow(s)) s[, `:=`(aln = j, strand = p$str, tp = p$tp)]
        s }))
      sg[, `:=`(qx0 = sqs + q$from, qx1 = sqe + q$from, tx0 = sts + t$from, tx1 = ste + t$from, identity = round(100 * match / (match + mism), 1))]
      sg[, covers := mapply(function(a, b) parts_of(segs, a, b), qx0, qx1)]
      setorder(sg, qx0)
      cat(sprintf("== alignment of %s's start to the DNA upstream of %s (minimap2 -x asm20): %d alignment(s), cut at gaps >= %d bp into %d segment(s)\n",
          ref$lab, tc$lab, nrow(paf), MINGAP, nrow(sg)))
      print(sg[, .(alignment = aln, `bp on complete copy (from its start codon)` = sprintf("%d..%d", as.integer(qx0), as.integer(qx1)),
                   `bp on incomplete copy (from its exon 2)` = sprintf("%d..%d", as.integer(tx0), as.integer(tx1)),
                   strand, type = fifelse(tp == "P", "primary", "secondary"), `identity %` = identity, covers)], row.names = FALSE)
      BL[[length(BL) + 1]] <- cbind(data.table(species = SP, copy = tc$lab, complete_copy = ref$lab), sg[, .(aln, qx0, qx1, tx0, tx1, strand, tp, identity, covers)])
    } else cat(sprintf("== alignment of %s's start to the DNA upstream of %s: no block\n", ref$lab, tc$lab))
    # where the complete copy's start codon and its exon 2 land (+ strand blocks, primary first)
    land <- function(xq) { if (!nrow(paf)) return(NA_real_); qo <- xq - q$from; cv <- paf[str == "+" & qs <= qo & qe > qo]
      if (!nrow(cv)) return(NA_real_); b <- cv[order(tp != "P", -nm)][1]; z <- cs_map(b$cs, b$qs, b$ts, qo); if (is.na(z)) NA_real_ else z + t$from }
    e2_x <- if (nrow(ref$ex) > 1) land(ref$ex$x0[2]) else NA_real_
    # the first exon itself, base by base
    pr <- substr(q$seq, ex1$x0 - q$from + 1, ex1$x1 - q$from)
    mmf <- ham(t$seq, pr); mmr <- ham(t$seq, rc(pr))
    comp <- table(factor(strsplit(t$seq, "")[[1]], levels = c("A", "C", "G", "T"))); pm <- sum((comp / sum(comp))^2)
    npos <- 2 * length(mmf)
    i1 <- which(mmf <= MAXMM); i2 <- which(mmr <= MAXMM)
    h <- data.table(strand = c(rep("same", length(i1)), rep("opposite", length(i2))), i = c(i1, i2), mm = c(mmf[i1], mmr[i2]))
    sub_t <- function(x0, x1) { a0 <- x0 - t$from + 1; a1 <- x1 - t$from; if (a0 < 1 || a1 > nchar(t$seq)) NA_character_ else substr(t$seq, a0, a1) }
    if (nrow(h)) {
      h[, `:=`(x = i - 1 + t$from, E = signif(npos * pbinom(mm, nchar(pr), 1 - pm), 2))]
      setorder(h, mm, -strand); h <- head(h, 5)
      h[, `:=`(sequence = mapply(function(x, s) { z <- sub_t(x, x + nchar(pr)); if (s == "same") z else rc(z) }, x, strand),
               genomic = sapply(x, function(z) tog(tc, z)))]
    } else h <- data.table(strand = character(), mm = integer(), x = numeric(), E = numeric(), sequence = character(), genomic = numeric())
    cat(sprintf("== %s's first exon (%d bp, %s) searched base by base over the same DNA, both strands: %d match(es) with <= %d mismatches\n",
        ref$lab, nchar(pr), pr, nrow(h), MAXMM))
    if (nrow(h)) print(h[, .(strand, mismatches = mm, `x on incomplete copy` = as.integer(x), genomic = sprintf("%s:%d", tc$contig, as.integer(genomic)),
                             E = formatC(E, format = "e", digits = 1), sequence)], row.names = FALSE)
    HI[[length(HI) + 1]] <- cbind(data.table(species = SP, copy = tc$lab), h[, .(strand, mismatches = mm, x, genomic, E, sequence)])
    aln_x <- land(ex1$x0)                                # where the complete copy's start codon aligns
    hs <- h[strand == "same" & E < 0.01]                 # best same-strand match; if several, the one at the aligned start codon
    hb <- if (nrow(hs) && !is.na(aln_x) && any(abs(hs$x - aln_x) <= 50)) hs[abs(x - aln_x) <= 50][1] else h[strand == "same"][1]
    probe_x <- if (nrow(h) && !is.na(hb$x) && hb$E < 0.01) hb$x else NA_real_
    src <- NA_character_; anc <- NA_real_
    if (!is.na(probe_x) && !is.na(aln_x)) { anc <- probe_x
      src <- if (abs(probe_x - aln_x) <= 50) "first-exon match and alignment agree" else
             sprintf("first-exon match (the alignment puts the start codon %d bp away)", as.integer(aln_x - probe_x))
    } else if (!is.na(probe_x)) { anc <- probe_x; src <- "first-exon match only (no alignment block covers it)"
    } else if (!is.na(aln_x)) { anc <- aln_x; src <- "alignment only (the first exon itself does not match)" }
    s_found <- if (is.na(anc)) NA_character_ else sub_t(anc, anc + L1)
    row <- data.table(species = SP, copy = tc$lab, contig = tc$contig, strand = tc$strand, aligned_from_aa = tc$first_aa,
      old_x0 = tc$g0, complete_copy = ref$lab, complete_contig = ref$contig, complete_strand = ref$strand, complete_x0 = ref$g0,
      first_exon_bp = L1, found = !is.na(anc), source = src, start_x = anc, start_genomic = if (is.na(anc)) NA_real_ else tog(tc, anc),
      exon2_lands_at_x = e2_x, alignment_start_x = aln_x, match_mismatches = if (nrow(h)) hb$mm else NA_integer_,
      match_E = if (nrow(h)) hb$E else NA_real_,
      start_codon = if (is.na(anc)) NA_character_ else sub_t(anc, anc + 3), donor = if (is.na(anc)) NA_character_ else sub_t(anc + L1, anc + L1 + 2),
      peptide_complete = translate(pr), peptide_found = translate(s_found),
      searched_from_x = t$from, genome_blocks_first_exon_or_upstream = NA_integer_,
      CpGs_in_window = NA_integer_, CpGs_with_calls = NA_integer_, median_reads = NA_real_, need_rescan = FALSE, bed = basename(bed_name(D)))
    if (!is.na(anc)) {
      afr <- list(contig = tc$contig, strand = tc$strand, g0 = tog(tc, anc))
      w <- gseq(afr, -2000, 500); m <- gregexpr("CG", w$seq)[[1]]; nref <- if (m[1] == -1) 0L else length(m)
      rr <- grange(afr, w$from, w$to)
      bb <- if (is.null(bed)) data.table() else bed[contig == afr$contig & pos >= rr[1] & pos <= rr[2]]
      row[, `:=`(CpGs_in_window = nref, CpGs_with_calls = nrow(bb), median_reads = if (nrow(bb)) median(bb$cov) else 0)]
      row[, need_rescan := CpGs_with_calls < 0.8 * CpGs_in_window | median_reads < 5]
      if (row$need_rescan) RS[[length(RS) + 1]] <- data.table(contig = afr$contig, lo = max(1, afr$g0 - FL), hi = min(clen(afr$contig), afr$g0 + FL))
      cat(sprintf("=> first exon of %s: %s:%d (%s strand), x = %d, i.e. %s upstream of its exon 2 (%s)\n",
          tc$lab, tc$contig, as.integer(afr$g0), tc$strand, as.integer(anc), kb(-anc), src))
      cat(sprintf("   start codon %s, splice donor after the exon %s; %s exon 1 encodes %s, %s's found exon %s; %s's exon 2 lands at x = %s\n",
          row$start_codon, row$donor, ref$lab, row$peptide_complete, tc$lab, row$peptide_found, ref$lab,
          if (is.na(e2_x)) "- (no block covers it)" else as.integer(e2_x)))
      cat(sprintf("   methylation calls already collected, -2 kb to +0.5 kb of it: %d of %d CpGs, median %s reads -> %s\n",
          row$CpGs_with_calls, row$CpGs_in_window, row$median_reads, if (row$need_rescan) "too thin: the shell step scans the BAMs again" else "enough"))
    } else {
      cat(sprintf("=> no start found for %s in the %s upstream of its exon 2; aligning %s's start to the whole genome (a few minutes)\n",
          tc$lab, kb(-t$from), ref$lab))
      wg <- paf_read(system2(MM2, c("-cx", "asm20", "--cs", "-N", "50", "-p", "0.05", "-s", "100", "-t", TH, shQuote(G), shQuote(qf)), stdout = TRUE, stderr = FALSE))
      if (nrow(wg)) {
        wg[, `:=`(qx0 = qs + q$from, qx1 = qe + q$from, identity = sapply(cs, cs_ident))][, covers := mapply(function(a, b) parts_of(segs, a, b), qx0, qx1)]
        wg <- wg[!(tname == ref$contig & te >= rq[1] & ts <= rq[2])]          # not the complete copy itself
      }
      fe <- if (nrow(wg)) wg[qx0 < ex1$x1] else wg                         # blocks reaching the first exon or the DNA upstream of it
      row[, genome_blocks_first_exon_or_upstream := nrow(fe)]
      if (nrow(wg)) { setorder(wg, -nm)
        cat(sprintf("   %d block(s) elsewhere in the genome (the complete copy itself left out), best 15:\n", nrow(wg)))
        print(head(wg, 15)[, .(where = sprintf("%s:%d-%d", tname, as.integer(ts + 1), as.integer(te)), strand = str, `identity %` = identity, covers)], row.names = FALSE)
        fwrite(cbind(data.table(species = SP, copy = tc$lab), wg[, .(tname, ts = ts + 1, te, strand = str, identity, covers)]),
               file.path(O, paste0("genome_blocks_", tc$lab, ".tsv")), sep = "\t")
      } else cat("   no block anywhere else in the genome\n")
      cat(sprintf("   => %s\n", if (nrow(fe)) sprintf("%d of them reach %s's first exon or the DNA upstream of it: see genome_blocks_%s.tsv", nrow(fe), ref$lab, tc$lab) else
          sprintf("nothing anywhere in the genome aligns to %s's first exon or the DNA upstream of it: %s has no recognisable first exon or promoter left (deleted or replaced)", ref$lab, tc$lab)))
    }
    AN[[length(AN) + 1]] <- row
  }
  an <- if (length(AN)) rbindlist(AN) else data.table(species = character(), copy = character(), found = logical(), need_rescan = logical())
  fwrite(an, file.path(O, "anchor.tsv"), sep = "\t")
  fwrite(if (length(BL)) rbindlist(BL) else data.table(species = character()), file.path(O, "blocks.tsv"), sep = "\t")
  fwrite(if (length(HI)) rbindlist(HI) else data.table(species = character()), file.path(O, "exon1_hits.tsv"), sep = "\t")
  fwrite(if (length(RS)) rbindlist(RS) else data.table(contig = character(), lo = integer(), hi = integer()), file.path(O, "rescan_region.tsv"), sep = "\t")
  bams <- if (file.exists(file.path(D, "bams.tsv"))) fread(file.path(D, "bams.tsv")) else data.table(n = integer(), bam = character())
  fwrite(bams[n %in% USE, .(n, bam)], file.path(O, "usable_bams.tsv"), sep = "\t")
  quit(save = "no")
}

# ---------------------------------------------------------------------------------------------------------- report
an <- fread(file.path(O, "anchor.tsv"))
if (!nrow(an)) { cat("   nothing to report: every TPXL6 copy is aligned from amino acid 1\n"); quit(save = "no") }
bed_main <- bed_read(bed_name(D)); bed_new <- bed_read(bed_name(O))
bl <- if (file.exists(file.path(O, "blocks.tsv")) && length(readLines(file.path(O, "blocks.tsv"), n = 2)) > 1) fread(file.path(O, "blocks.tsv")) else data.table()
bams <- if (file.exists(file.path(D, "bams.tsv"))) fread(file.path(D, "bams.tsv")) else data.table(n = integer(), bam = character())
ctl_row <- L[species == SP & query == Q5][order(-score)][1]
CTL <- if (nrow(ctl_row) && !is.na(ctl_row$contig)) copy_of(ctl_row, Q5) else NULL

# rows of the methylation panel: the complete copy, each incomplete copy from its found start, the related gene
R <- list(list(lab = paste("TPXL6", REF$lab, "(complete copy)"), fr = REF[c("contig", "strand", "g0")], ex = REF$ex[, .(x0, x1)],
               red = NULL, bed = bed_main, what = sprintf("its first exon (miniprot), %s:%d", REF$contig, as.integer(REF$g0))))
for (i in seq_len(nrow(an))) {
  r <- an[i]; tc <- CC[[which(sapply(CC, `[[`, "contig") == r$contig & sapply(CC, `[[`, "g0") == r$old_x0)[1]]]
  if (!isTRUE(r$found)) { R[[length(R) + 1]] <- list(lab = paste("TPXL6", r$copy, "from exon 2 (no first exon found)"), fr = tc[c("contig", "strand", "g0")], ex = tc$ex[, .(x0, x1)],
      red = NULL, bed = bed_main, what = "exon 2 (no first exon found)", nf = TRUE); next }
  R[[length(R) + 1]] <- list(lab = sprintf("TPXL6 %s, from the first exon found %s upstream of exon 2", r$copy, kb(-r$start_x)),
    fr = list(contig = r$contig, strand = r$strand, g0 = r$start_genomic), ex = tc$ex[, .(x0 = x0 - r$start_x, x1 = x1 - r$start_x)],
    red = data.table(x0 = 0, x1 = r$first_exon_bp), bed = if (is.null(bed_new)) bed_main else bed_new,
    what = sprintf("the found first exon, %s:%d (%s)", r$contig, as.integer(r$start_genomic), r$source))
}
if (!is.null(CTL)) R[[length(R) + 1]] <- list(lab = paste0("related gene, ", CTL$lab, " (expressed)"), fr = CTL[c("contig", "strand", "g0")],
  ex = CTL$ex[, .(x0, x1)], red = NULL, bed = bed_main, what = sprintf("its first exon (miniprot), %s:%d", CTL$contig, as.integer(CTL$g0)))
rng <- rbindlist(lapply(R, function(z) { gr <- grange(z$fr, -6000, 6000); data.table(contig = z$fr$contig, lo = gr[1], hi = gr[2]) }))
rp <- carp_get(rng); if (!is.null(rp)) rp <- rp[!is.na(grp)]
out <- list()
for (k in seq_along(R)) {
  z <- R[[k]]; fr <- z$fr; gr <- grange(fr, -6000, 6000)
  rr <- if (is.null(rp)) data.table(grp = character(), x0 = numeric(), x1 = numeric()) else
        rp[contig == fr$contig & end >= gr[1] & start <= gr[2], .(grp, x0 = pmin(tox(fr, start), tox(fr, end)), x1 = pmax(tox(fr, start), tox(fr, end)) + 1)]
  te <- rr[grp != "rep"]
  cg <- if (is.null(z$bed)) data.table(x = numeric(), mod = numeric(), cov = numeric()) else
        z$bed[contig == fr$contig & pos >= gr[1] & pos <= gr[2], .(x = tox(fr, pos), mod, cov)][order(x)]
  cg[, in_te := inside(x, te$x0, te$x1)]
  if (nrow(cg)) cg[, smooth := sapply(x, function(v) mean(mod[abs(x - v) <= SMOOTH]))]
  res <- data.table(species = SP, region = z$lab, start_used = z$what)
  for (w in names(WIN)) { v <- cg[!in_te & x >= WIN[[w]][1] & x < WIN[[w]][2]]
    res[[paste0("meth_", w)]] <- if (nrow(v)) round(mean(v$mod), 1) else NA_real_; res[[paste0("n_", w)]] <- nrow(v)
    res[[paste0("reads_", w)]] <- if (nrow(v)) median(v$cov) else NA_real_ }
  out[[k]] <- res
  R[[k]]$rr <- rr[x1 > X0 & x0 < X1]; R[[k]]$cg <- cg[x >= X0 & x <= X1]
}
tab <- rbindlist(out, fill = TRUE)
fwrite(tab, paste0(OUTP, ".tsv"), sep = "\t")

# ---- figure: A how the complete copy's start aligns to the incomplete copy; B methylation around each start
r1 <- an[1]; tc1 <- CC[[which(sapply(CC, `[[`, "contig") == r1$contig & sapply(CC, `[[`, "g0") == r1$old_x0)[1]]]
found1 <- isTRUE(r1$found)
sR <- if (found1) 0 else REF$ex$x0[2]; sT <- if (found1) r1$start_x else 0        # plotted x = gene x - shift
qa <- REF$ex$x0[1] - UP; qb <- if (nrow(REF$ex) > 1) REF$ex$x1[2] + 200 else REF$ex$x1[1] + 2000
ta <- if (found1) r1$start_x - UP else -UP; tb <- tc1$ex$x1[1] + 200
b1 <- if (nrow(bl)) bl[copy == r1$copy] else bl
xa0 <- min(qa - sR, ta - sT); xa1 <- max(qb - sR, tb - sT)
far <- if (nrow(b1)) b1[tx1 - sT < xa0 - 2000 | tx0 - sT > xa1 + 2000] else b1
b1 <- if (nrow(b1)) b1[!(tx1 - sT < xa0 - 2000 | tx0 - sT > xa1 + 2000)] else b1
if (nrow(b1)) { xa0 <- min(xa0, b1$tx0 - sT, b1$qx0 - sR); xa1 <- max(xa1, b1$tx1 - sT, b1$qx1 - sR) }
ga <- rbind(data.table(contig = REF$contig, lo = grange(REF, qa, qb)[1], hi = grange(REF, qa, qb)[2]),
            data.table(contig = tc1$contig, lo = grange(tc1, xa0 + sT, xa1 + sT)[1], hi = grange(tc1, xa0 + sT, xa1 + sT)[2]))
rpa <- carp_get(ga); if (!is.null(rpa)) rpa <- rpa[!is.na(grp)]
tes_of <- function(fr, s) if (is.null(rpa)) data.table() else
  rpa[contig == fr$contig, .(grp, x0 = pmin(tox(fr, start), tox(fr, end)) - s, x1 = pmax(tox(fr, start), tox(fr, end)) + 1 - s)][x1 > xa0 & x0 < xa1]
nB <- length(R); HA <- 2.9; HB <- 1.45 + 0.95 * nB
draw <- function() {
  layout(matrix(1:2, ncol = 1), heights = c(HA, HB))
  # A
  par(mai = c(0.7, 2.3, 1.05, 0.3), xpd = NA, family = "sans")
  plot(NA, xlim = c(xa0, xa1), ylim = c(2.5, 0.6), xaxs = "i", axes = FALSE, xlab = "", ylab = "")
  span <- xa1 - xa0; minw <- 0.0035 * span
  at <- pretty(c(xa0, xa1) / 1000); at <- at[at * 1000 >= xa0 & at * 1000 <= xa1]
  axis(1, at = at * 1000, labels = at, cex.axis = 0.72, col = "grey55", col.axis = INK2)
  mtext(if (found1) "kb from the start codon (top: its first exon; bottom: the first exon found here)" else
        "kb from exon 2 (no first exon found on the bottom copy)", side = 1, line = 2.1, cex = 0.68, col = INK2)
  yT <- 1; yB <- 2; hh <- 0.09
  if (nrow(b1)) for (j in seq_len(nrow(b1))) { b <- b1[j]
    xs <- if (b$strand == "+") c(b$qx0 - sR, b$qx1 - sR, b$tx1 - sT, b$tx0 - sT) else c(b$qx0 - sR, b$qx1 - sR, b$tx0 - sT, b$tx1 - sT)
    polygon(xs, c(yT + hh, yT + hh, yB - hh, yB - hh), col = alpha(RIB, if (b$tp == "P") 0.35 else 0.15), border = NA)
    if ((b$qx1 - b$qx0) > 0.05 * span) text(mean(xs), (yT + yB) / 2, sprintf("%.0f%%", b$identity), cex = 0.6, col = INK2) }
  if (found1 && grepl("first-exon match", r1$source)) {
    polygon(c(0, max(r1$first_exon_bp, minw), max(r1$first_exon_bp, minw), 0), c(yT + hh, yT + hh, yB - hh, yB - hh), col = alpha(RED, 0.75), border = NA)
    text(0, yT - 0.2, sprintf("first exon matched base by base: %d mismatch%s in %d bp", r1$match_mismatches,
         if (r1$match_mismatches == 1) "" else "es", r1$first_exon_bp), adj = c(0, 0.5), cex = 0.62, col = RED) }
  track <- function(y, fr, s, ex, a0, a1, red = NULL) {
    a0 <- max(a0, xa0); a1 <- min(a1, xa1)
    segments(a0, y, a1, y, col = "#d6d5d0", lwd = 1.2)
    tt <- tes_of(fr, s); if (nrow(tt)) tt <- tt[x1 > a0 & x0 < a1]
    if (nrow(tt)) rect(pmax(tt$x0, a0), y - hh, pmin(tt$x1, a1), y + hh, col = sapply(GCOL[tt$grp], alpha, 0.45), border = NA)
    e <- ex[x1 - s > a0 & x0 - s < a1]
    if (nrow(e)) rect(e$x0 - s, y - 0.06, pmax(e$x1 - s, e$x0 - s + minw), y + 0.06, col = INK, border = NA)
    if (!is.null(red)) rect(red$x0, y - 0.06, pmax(red$x1, red$x0 + minw), y + 0.06, col = RED, border = NA)
  }
  track(yT, REF, sR, REF$ex, qa - sR, qb - sR)
  track(yB, tc1, sT, tc1$ex, xa0, tb - sT, if (found1) data.table(x0 = 0, x1 = r1$first_exon_bp) else NULL)
  text(xa0 - 0.02 * span, yT, sprintf("%s TPXL6 %s\ncomplete copy (%s strand)", nice(SP), REF$lab, REF$strand), adj = 1, cex = 0.72, col = INK)
  text(xa0 - 0.02 * span, yB, sprintf("%s TPXL6 %s\nminiprot from amino acid %d (%s strand)", nice(SP), r1$copy, r1$aligned_from_aa, r1$strand), adj = 1, cex = 0.72, col = INK)
  if (found1) text(0, yB + 0.17, "found\nfirst exon", col = RED, cex = 0.6, adj = c(0.5, 1))
  text(tc1$ex$x0[1] - sT, yB + 0.17, "exon 2\n(miniprot)", col = INK, cex = 0.6, adj = c(0.5, 1))
  mtext(sprintf("A  Where the start of %s aligns on %s", REF$lab, r1$copy), side = 3, line = 3.75, adj = 0, at = xa0 - 0.02 * span - 0.0, cex = 1.0, font = 2, col = INK)
  mtext(sprintf("Grey: segments aligned by minimap2 (-x asm20; repeats masked on the top copy; cut at gaps >= %d bp), %% identity.%s", MINGAP,
                if (nrow(far)) sprintf(" %d more segment(s) further away, not drawn (blocks.tsv).", nrow(far)) else ""),
        side = 3, line = 2.75, adj = 0, cex = 0.64, col = INK2)
  mtext(if (found1 && grepl("first-exon match", r1$source)) sprintf("Red: the top copy's first exon matched base by base on the bottom copy (E = %s).", format(r1$match_E, digits = 2, scientific = TRUE))
        else if (found1) "Red: where the top copy's start codon aligns (its first exon itself does not match base by base)."
        else sprintf("No match for the top copy's first exon (<= %d mismatches) in the %s searched upstream of the bottom copy%s.", MAXMM, kb(-r1$searched_from_x),
                     if (isTRUE(r1$genome_blocks_first_exon_or_upstream == 0)) ", nor anywhere else in the genome" else ""),
        side = 3, line = 2.0, adj = 0, cex = 0.64, col = if (found1) RED else INK2)
  mtext("Black: coding exons (miniprot). Shaded: CARP repeats (blue LTR, pink other retrotransposons, yellow DNA transposons, grey other).",
        side = 3, line = 1.25, adj = 0, cex = 0.64, col = INK2)
  if (nzchar(Sys.getenv("FIG_BANNER"))) mtext(Sys.getenv("FIG_BANNER"), side = 3, line = 3.75, adj = 1, cex = 0.95, font = 2, col = RED)
  # B
  par(mai = c(0.7, 2.3, 0.85, 0.3), xpd = NA)
  plot(NA, xlim = c(X0, X1), ylim = c(nB + 0.6, 0.45), xaxs = "i", axes = FALSE, xlab = "", ylab = "")
  axis(1, at = seq(X0, X1, 500), labels = seq(X0, X1, 500) / 1000, cex.axis = 0.72, col = "grey55", col.axis = INK2)
  mtext("kb from each gene's start codon (x = 0)", side = 1, line = 2.2, cex = 0.7, col = INK2)
  for (k in seq_len(nB)) {
    z <- R[[k]]; yt <- k - 0.40; yb <- k + 0.12; yg <- k + 0.30; ytr <- function(m) yb - m / 100 * (yb - yt)
    text(X0 - 0.03 * (X1 - X0), k - 0.2, paste(strwrap(z$lab, 34), collapse = "\n"), adj = 1, cex = 0.7, col = if (is.null(z$red)) INK else RED)
    if (nrow(z$rr)) rect(pmax(z$rr$x0, X0), yt, pmin(z$rr$x1, X1), yb, col = sapply(GCOL[z$rr$grp], alpha, 0.28), border = NA)
    segments(X0, c(yt, ytr(50), yb), X1, c(yt, ytr(50), yb), col = c("grey88", "grey80", "grey70"), lty = c(1, 3, 1), lwd = 0.7)
    text(X0 - 0.005 * (X1 - X0), c(yt, yb), c("100%", "0"), adj = c(1, 0.5), cex = 0.5, col = MUTED)
    if (nrow(z$cg)) { points(z$cg$x, ytr(z$cg$mod), pch = 16, cex = 0.35, col = alpha(INK, 0.5)); lines(z$cg$x, ytr(z$cg$smooth), col = INK, lwd = 1.2) }
    else text((X0 + X1) / 2, (yt + yb) / 2, "no methylation calls here", cex = 0.7, font = 3, col = MUTED)
    segments(X0, yg, X1, yg, col = "#d6d5d0"); e <- z$ex[x1 > X0 & x0 < X1]
    if (nrow(e)) rect(pmax(e$x0, X0), yg - 0.07, pmin(e$x1, X1), yg + 0.07, col = INK, border = NA)
    if (!is.null(z$red)) rect(z$red$x0, yg - 0.07, pmax(z$red$x1, z$red$x0 + 15), yg + 0.07, col = RED, border = NA)
    segments(0, yt, 0, yg, lty = 3, col = MUTED)
  }
  mtext("B  CpG methylation around each start (same BAM for all rows)", side = 3, line = 2.9, adj = 0, at = X0 - 0.03 * (X1 - X0) - 0.0, cex = 1.0, font = 2, col = INK)
  mtext(sprintf("Dots: CpGs (%% of reads methylated, >= 3 reads); line: mean within +-%d bp. Black: coding exons (miniprot); red: the first exon found in panel A; dotted: x = 0.", SMOOTH),
        side = 3, line = 1.7, adj = 0, cex = 0.64, col = INK2)
}
H_IN <- HA + HB
dev_png <- function(f, w, h) if (capabilities("cairo")) png(f, width = w, height = h, units = "in", res = 220, type = "cairo") else png(f, width = w, height = h, units = "in", res = 220)
dev_png(paste0(OUTP, ".png"), 10.5, H_IN); draw(); invisible(dev.off())
pdf(paste0(OUTP, ".pdf"), 10.5, H_IN); draw(); invisible(dev.off())

ubam <- bams[n %in% USE]$bam
cat(sprintf("== start used for each row (BAM%s: %s)\n", if (length(ubam) > 1) "s" else "", paste(basename(ubam), collapse = ", ")))
for (z in R) cat(sprintf("   %-70s %s\n", z$lab, z$what))
compact <- function(d, lab) { p <- data.table(region = lab)
  for (w in names(WIN)) { m <- d[[paste0("meth_", w)]]; n <- d[[paste0("n_", w)]]
    p[[paste(w, "bp")]] <- ifelse(is.na(m), "-", sprintf("%.1f (%d)", m, as.integer(n))) }
  rd <- grep("^reads_", names(d), value = TRUE)
  if (length(rd)) p[["reads per CpG"]] <- apply(d[, ..rd], 1, function(v) if (all(is.na(v))) "-" else
                                                paste(unique(range(round(v), na.rm = TRUE)), collapse = "-"))
  p }
cat("\n== mean CpG methylation (%) by window around each start, CpGs inside TEs left out; (number of CpGs); median reads per CpG\n")
print(compact(tab, tab$region), row.names = FALSE)
tss <- file.path(H, "figures/tpxl6_methyl_tss_check.tsv")
if (file.exists(tss)) { ts <- fread(tss); oth <- ts[species != SP]
  if (nrow(oth)) { cat("\n== for reference, the other species from tpxl6_methyl_tss_check.tsv (their own BAMs)\n")
    print(compact(oth, paste(nice(oth$species), oth$region)), row.names = FALSE) } }
cat(sprintf("wrote %s.png / .pdf / .tsv\n", OUTP))
