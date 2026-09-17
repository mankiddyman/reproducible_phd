# ---------------------------------------------------------------------------
# DR18b -- identify GENESPACE's plotted gene subset, then map bp -> plot x.
# Read-only on inputs. Writes DR/out/DR18b_gene_x.csv and DR18b_chr_map.csv
# base R + data.table only (genespace module = R 4.2.0, no dplyr)
# ---------------------------------------------------------------------------
suppressMessages(library(data.table))
options(width = 210)

RDA      <- "../genespace/riparian/Nepenthes_gracilis_geneOrder_rSourceData.rda"
COMBBED  <- "../genespace/results/combBed.txt"
OUT_GENE <- "DR/out/DR18b_gene_x.csv"
OUT_CHR  <- "DR/out/DR18b_chr_map.csv"
stopifnot(file.exists(RDA), file.exists(COMBBED))

e <- new.env(); load(RDA, envir = e)
P     <- e$srcd$ggplotObj
blk   <- as.data.table(e$srcd$sourceData$blocks)
chrom <- as.data.table(e$srcd$sourceData$chromosomes)
br    <- as.data.table(P$layers[[1]]$data)
cb    <- fread(COMBBED)
CH    <- chrom[, .(genome, chr, len = length, x1, x2)]

# --- 1. FLAG COLUMN STRUCTURE ----------------------------------------------
cat("\n================ 1. FLAG COLUMNS ================\n")
for (nm in c("isArrayRep","noAnchor")) { cat("--", nm, ":\n"); print(table(cb[[nm]], useNA = "ifany")) }
cat("arrayID NA or empty :", sum(is.na(cb$arrayID) | cb$arrayID == ""), "of", nrow(cb), "\n")
cat("og NA               :", sum(is.na(cb$og)), "\n")
cat("globOG NA or empty  :", sum(is.na(cb$globOG)  | cb$globOG  == ""), "\n")
cat("globHOG NA or empty :", sum(is.na(cb$globHOG) | cb$globHOG == ""), "\n")
cat("total combBed deficit vs sum(length):",
    nrow(cb[CH, on = .(genome, chr), nomatch = 0]) - sum(CH$len), "\n")

# --- 2. ANCHORS: genes whose plot ord is known -----------------------------
anch <- unique(rbind(
  blk[, .(genome = genome1, chr = chr1, gene = firstGene1, pord = startOrd1)],
  blk[, .(genome = genome1, chr = chr1, gene = lastGene1,  pord = endOrd1)],
  blk[, .(genome = genome2, chr = chr2, gene = firstGene2, pord = startOrd2)],
  blk[, .(genome = genome2, chr = chr2, gene = lastGene2,  pord = endOrd2)]))
anch <- merge(anch, CH[, .(genome, chr)], by = c("genome","chr"))
cat("\nanchor genes with known plot ord:", nrow(anch),
    "| chromosomes covered:", uniqueN(anch[, .(genome, chr)]), "\n")
cat("anchors per chromosome:\n"); print(summary(anch[, .N, by = .(genome, chr)]$N))
dupg <- anch[, .N, by = .(genome, chr, gene)][N > 1]
cat("genes appearing with two different plot ords:", nrow(anch[, uniqueN(pord), by=.(genome,chr,gene)][V1>1]), "\n")

# --- 3. CANDIDATE SUBSETS, TESTED BY RANK ----------------------------------
cat("\n================ 3. SUBSET SEARCH ================\n")
keep <- list(
  all               = rep(TRUE, nrow(cb)),
  arrayRepStrict    = cb$isArrayRep %in% TRUE,
  notArrayDup       = !(cb$isArrayRep %in% FALSE),
  anchorable        = !(cb$noAnchor %in% TRUE),
  hasOG             = !is.na(cb$og),
  hasGlobOG         = !is.na(cb$globOG)  & cb$globOG  != "",
  hasGlobHOG        = !is.na(cb$globHOG) & cb$globHOG != "",
  notArrayDup_anch  = !(cb$isArrayRep %in% FALSE) & !(cb$noAnchor %in% TRUE),
  notArrayDup_hasOG = !(cb$isArrayRep %in% FALSE) & !is.na(cb$og))

rankset <- function(mask) {
  gg <- cb[mask, .(genome, chr, ofID, start, end)]
  gg <- merge(gg, CH, by = c("genome","chr"))
  setorder(gg, genome, chr, start)
  gg[, r  := seq_len(.N), by = .(genome, chr)]
  gg[, NN := .N,          by = .(genome, chr)]
  gg }

score <- function(gg) {
  a <- merge(anch, gg[, .(genome, chr, gene = ofID, r, NN, len)],
             by = c("genome","chr","gene"), all.x = TRUE)
  a[, `:=`(fw = r == pord, rv = (NN + 1L - r) == pord)]
  s <- a[, .(n = .N, unres = sum(is.na(r)),
             fw = sum(fw, na.rm = TRUE), rv = sum(rv, na.rm = TRUE),
             NN = NN[1], len = len[1]), by = .(genome, chr)]
  s[, mode := fifelse(unres == 0 & fw == n, "fwd",
              fifelse(unres == 0 & rv == n, "rev", NA_character_))]
  s }

best <- NA_character_; tab <- list()
for (nm in names(keep)) {
  s <- score(rankset(keep[[nm]])); tab[[nm]] <- s
  cat(sprintf("  %-18s N==length %2d/87 | anchors consistent %2d/87 (fwd %2d rev %2d) | unresolved anchors %5d\n",
      nm, sum(s$NN == s$len, na.rm = TRUE), sum(!is.na(s$mode)),
      sum(s$mode == "fwd", na.rm = TRUE), sum(s$mode == "rev", na.rm = TRUE), sum(s$unres)))
  if (is.na(best) && sum(!is.na(s$mode)) == nrow(s)) best <- nm }

if (is.na(best)) {
  cat("\n================ FORENSIC: WHICH GENES ARE DROPPED? ================\n")
  gg <- rankset(keep$all)
  a  <- merge(anch, gg[, .(genome, chr, gene = ofID, r, NN, len)],
              by = c("genome","chr","gene"), all.x = TRUE)
  a[, d := r - pord]
  setorder(a, genome, chr, r)
  cat("drift d = rank_all - plot_ord, by chromosome (monotone rise => genes dropped along chr):\n")
  print(a[, .(minD = min(d), maxD = max(d), deficit = NN[1] - len[1],
              spearman = suppressWarnings(cor(r, d, method = "spearman"))),
          by = .(genome, chr)][order(-deficit)][1:15])
  tg <- a[, .(deficit = NN[1] - len[1]), by = .(genome, chr)][deficit %in% 1:3][1:3]
  for (i in seq_len(nrow(tg))) {
    aa <- a[tg[i], on = .(genome, chr)][order(r)]
    j  <- which(diff(aa$d) != 0)
    cat("\n--", tg$genome[i], tg$chr[i], "-- deficit", tg$deficit[i],
        "-- d jumps between anchor ranks:", paste(aa$r[j], aa$r[j+1], sep="-", collapse=", "), "\n")
    if (length(j)) {
      w <- gg[.(tg$genome[i], tg$chr[i]), on = .(genome, chr)][r > aa$r[j[1]] & r < aa$r[j[1]+1]]
      cat("   candidate dropped genes in that window (", nrow(w), " genes):\n", sep="")
      print(merge(w[, .(ofID, r, start)], cb[, .(ofID, id, pepLen, arrayID, isArrayRep, noAnchor, og, globOG, globHOG)],
                  by = "ofID")[order(r)], nrows = 40) } }
  stop("no candidate subset reproduces GENESPACE's plot ord -- inspect the windows above")
}

cat("\n>>> WINNING SUBSET:", best, "\n")
s <- tab[[best]]
cat("chromosomes where N == length:", sum(s$NN == s$len), "/", nrow(s), "\n")
if (any(s$NN != s$len)) print(s[NN != len][order(genome, chr)])
cat("\nflip state:\n"); print(table(s$mode))
cat("flipped chromosomes:", sum(s$mode == "rev"), "\n")
if (any(s$mode == "rev")) print(s[mode == "rev", .(genome, chr, n, NN, len)][order(genome, chr)], nrows = 100)
cat("\nthe chromosome you flagged:\n"); print(s[genome=="Drosera_capensis" & chr=="chr15_collapsed"])

# --- 4. GENE -> PLOT X ------------------------------------------------------
cat("\n================ 4. GENE -> PLOT X ================\n")
g <- rankset(keep[[best]])
g <- merge(g, s[, .(genome, chr, flipped = mode == "rev")], by = c("genome","chr"))
g[, plot_ord := fifelse(flipped, NN + 1L - r, r)]
g[, x := x1 + plot_ord]
cat("every gene x inside its chromosome box? ", all(g$x >= g$x1 & g$x <= g$x2), "\n")
print(g[, .(lo = min(x), hi = max(x), x1 = x1[1], x2 = x2[1], flipped = flipped[1]),
        by = .(genome, chr)][order(genome, chr)][1:6])

# --- 5. VALIDATION AGAINST THE ACTUAL POLYGONS -----------------------------
cat("\n================ 5. VALIDATION ================\n")
setorder(g, genome, chr, start)
K <- g[!duplicated(g[, .(genome, chr, start)])]
mb2x <- function(gn, ch, bp) {
  d <- K[.(gn, ch), on = .(genome, chr)]
  if (!nrow(d)) return(rep(NA_real_, length(bp)))
  approx(c(0, d$start, max(d$end) + 1),
         c(if (d$flipped[1]) d$x2[1] else d$x1[1], d$x,
           if (d$flipped[1]) d$x1[1] else d$x2[1]), xout = bp, rule = 2)$y }
sides <- rbind(
  blk[, .(blkID, yrow = as.numeric(y1), genome = genome1, chr = chr1, bpa = startBp1, bpb = endBp1)],
  blk[, .(blkID, yrow = as.numeric(y2), genome = genome2, chr = chr2, bpa = startBp2, bpb = endBp2)])
ends  <- br[, .(xlo = min(x), xhi = max(x)), by = .(blkID, yrow = round(y, 6))]
sides <- merge(sides, ends, by = c("blkID","yrow"), all.x = TRUE)
cat("block sides matched to a polygon attachment:", sum(!is.na(sides$xlo)), "/", nrow(sides), "\n")
sides[, pa := mb2x(genome[1], chr[1], bpa), by = .(genome, chr)]
sides[, pb := mb2x(genome[1], chr[1], bpb), by = .(genome, chr)]
sides[, resid := pmax(abs(pmin(pa,pb) - xlo), abs(pmax(pa,pb) - xhi))]
cat("\nresidual |predicted x - GENESPACE polygon x|, in genes:\n"); print(summary(sides$resid))
cat("\nworst by genome:\n")
print(sides[, .(max_resid = max(resid, na.rm=TRUE), n = .N), by = genome][order(-max_resid)])
cat("\n10 worst block sides:\n")
print(sides[order(-resid)][1:10, .(genome, chr, bpa, bpb, pa, pb, xlo, xhi, resid)])

# --- 6. WRITE + VERDICT ----------------------------------------------------
fwrite(g[, .(genome, chr, ofID, start, end, r, NN, len, flipped, plot_ord, x1, x2, x)], OUT_GENE)
fwrite(merge(s[, .(genome, chr, n_anchors = n, mode, nSubset = NN, nGenes = len)],
             chrom[, .(genome, chr, x1, x2)], by = c("genome","chr")), OUT_CHR)
mx <- max(sides$resid, na.rm = TRUE)
cat("\n################ VERDICT ################\n")
cat("subset:", best, "| flipped chromosomes:", sum(s$mode == "rev"), "/", nrow(s), "\n")
cat("max residual vs GENESPACE polygons:", signif(mx, 5), "genes\n")
cat("sides with residual > 2 genes:", sum(sides$resid > 2, na.rm = TRUE), "\n")
if (is.na(mx) || mx > 5) { cat("\n>>> FAILED: bp -> x does not reproduce GENESPACE's attachments.\n")
  stop("mapping rejected -- inspect the 10 worst sides above") }
cat("\n>>> bp -> x REPRODUCES GENESPACE's own attachments to within", signif(mx,3), "genes.\n")
cat("written:", OUT_GENE, "|", OUT_CHR, "\n")
