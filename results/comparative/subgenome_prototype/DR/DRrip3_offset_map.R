# ---------------------------------------------------------------------------
# DR18c -- plot_ord = combBed rank - per-chromosome offset. Build + validate.
# Read-only on inputs. Writes DR/out/DR18c_gene_x.csv, DR18c_chr_map.csv
# base R + data.table only (genespace module = R 4.2.0, no dplyr)
# ---------------------------------------------------------------------------
suppressMessages(library(data.table))
options(width = 210)

RDA      <- "../genespace/riparian/Nepenthes_gracilis_geneOrder_rSourceData.rda"
COMBBED  <- "../genespace/results/combBed.txt"
OUT_GENE <- "DR/out/DR18c_gene_x.csv"
OUT_CHR  <- "DR/out/DR18c_chr_map.csv"
stopifnot(file.exists(RDA), file.exists(COMBBED))

e <- new.env(); load(RDA, envir = e)
blk   <- as.data.table(e$srcd$sourceData$blocks)
chrom <- as.data.table(e$srcd$sourceData$chromosomes)
br    <- as.data.table(e$srcd$ggplotObj$layers[[1]]$data)
cb    <- fread(COMBBED)
CH    <- chrom[, .(genome, chr, len = length, x1, x2)]

# --- 1. RANK ----------------------------------------------------------------
cat("\n================ 1. RANK ================\n")
g <- merge(cb[, .(genome, chr, ofID, id, start, end, ord)], CH, by = c("genome","chr"))
setorder(g, genome, chr, ord)
g[, r := ord - min(ord) + 1L, by = .(genome, chr)]
cat("genes:", nrow(g), "| chromosomes:", uniqueN(g[, .(genome, chr)]), "\n")
cat("rank by ord == rank by start on every chromosome? ",
    all(g[, .(ok = identical(r, as.integer(frank(start, ties.method = 'first')))),
          by = .(genome, chr)]$ok), "\n")

# --- 2. OFFSET FROM ANCHORS -------------------------------------------------
cat("\n================ 2. OFFSET ================\n")
A <- unique(rbind(
  blk[, .(genome = genome1, chr = chr1, gene = firstGene1, pord = startOrd1)],
  blk[, .(genome = genome1, chr = chr1, gene = lastGene1,  pord = endOrd1)],
  blk[, .(genome = genome2, chr = chr2, gene = firstGene2, pord = startOrd2)],
  blk[, .(genome = genome2, chr = chr2, gene = lastGene2,  pord = endOrd2)]))
A <- merge(A, g[, .(genome, chr, gene = ofID, r)], by = c("genome","chr","gene"))
A[, d := r - pord]
setorder(A, genome, chr, r)
cat("anchors resolved:", nrow(A), "| chromosomes:", uniqueN(A[, .(genome, chr)]), "\n")

sm <- A[, .(nAnch = .N, minD = min(d), maxD = max(d), nD = uniqueN(d)), by = .(genome, chr)]
cat("\noffset spread per chromosome (maxD - minD):\n"); print(table(sm$maxD - sm$minD))
cat("\noffset magnitude:\n"); print(summary(sm$minD))
cat("chromosomes needing a step function (spread > 0):", sum(sm$maxD > sm$minD), "\n")
if (any(sm$maxD > sm$minD)) print(sm[maxD > minD][order(genome, chr)])
if (any(sm$maxD - sm$minD > 3)) { print(sm[maxD - minD > 3][order(-(maxD-minD))])
  stop("offset is not near-constant on some chromosome -- see above") }

setkey(A, genome, chr, r)
g[, off := { aa <- A[.(.BY$genome, .BY$chr)]
             aa$d[pmax(1L, findInterval(r, aa$r))] }, by = .(genome, chr)]
g[, plot_ord := r - off]
cat("\nplot_ord strictly increasing within every chromosome? ",
    all(g[, .(ok = !is.unsorted(plot_ord, strictly = TRUE)), by = .(genome, chr)]$ok), "\n")
cat("genes below the plotted range (plot_ord < 1):  ", sum(g$plot_ord < 1), "\n")
cat("genes above the plotted range (plot_ord > len):", sum(g$plot_ord > g$len), "\n")
cat("  (these are the leading/trailing genes GENESPACE dropped -- clamped to the box edge)\n")

# --- 3. GENE -> X -----------------------------------------------------------
g[, x := x1 + pmin(pmax(plot_ord, 1L), len)]
cat("\nevery x inside its chromosome box? ", all(g$x >= g$x1 & g$x <= g$x2), "\n")
print(g[, .(plotted = sum(plot_ord >= 1 & plot_ord <= len), len = len[1],
            off = off[1], lo = min(x), hi = max(x), x1 = x1[1], x2 = x2[1]),
        by = .(genome, chr)][order(genome, chr)][1:8])

# --- 4. VALIDATION AGAINST THE ACTUAL POLYGONS ------------------------------
cat("\n================ 4. VALIDATION ================\n")
setorder(g, genome, chr, start)
K <- g[!duplicated(g[, .(genome, chr, start)])]
setkey(K, genome, chr)
mb2x <- function(gn, ch, bp) { d <- K[.(gn, ch)]
  approx(c(0, d$start, max(d$end) + 1), c(d$x1[1], d$x, d$x2[1]), xout = bp, rule = 2)$y }
sides <- rbind(
  blk[, .(blkID, yrow = as.numeric(y1), genome = genome1, chr = chr1, bpa = startBp1, bpb = endBp1)],
  blk[, .(blkID, yrow = as.numeric(y2), genome = genome2, chr = chr2, bpa = startBp2, bpb = endBp2)])
sides <- merge(sides, br[, .(xlo = min(x), xhi = max(x)), by = .(blkID, yrow = round(y, 6))],
               by = c("blkID","yrow"), all.x = TRUE)
cat("block sides matched to a polygon attachment:", sum(!is.na(sides$xlo)), "/", nrow(sides), "\n")
sides[, pa := mb2x(genome[1], chr[1], bpa), by = .(genome, chr)]
sides[, pb := mb2x(genome[1], chr[1], bpb), by = .(genome, chr)]
sides[, resid := pmax(abs(pmin(pa, pb) - xlo), abs(pmax(pa, pb) - xhi))]
cat("\nresidual |predicted x - GENESPACE polygon x|, in genes:\n"); print(summary(sides$resid))
cat("\nworst by genome:\n")
print(sides[, .(max_resid = max(resid, na.rm = TRUE), n = .N), by = genome][order(-max_resid)])
cat("\n10 worst block sides:\n")
print(sides[order(-resid)][1:10, .(genome, chr, bpa, bpb, pa, pb, xlo, xhi, resid)])

# --- 5. WRITE + VERDICT -----------------------------------------------------
fwrite(g[, .(genome, chr, ofID, id, start, end, r, off, plot_ord, len, x1, x2, x)], OUT_GENE)
fwrite(merge(sm, chrom[, .(genome, chr, nGenes = length, x1, x2)], by = c("genome","chr")), OUT_CHR)
mx <- max(sides$resid, na.rm = TRUE)
cat("\n################ VERDICT ################\n")
cat("flipped chromosomes: 0 (offset is constant, a flip would not be)\n")
cat("max residual vs GENESPACE polygons:", signif(mx, 5), "genes\n")
cat("sides with residual > 1 gene:", sum(sides$resid > 1, na.rm = TRUE), "/", nrow(sides), "\n")
if (is.na(mx) || mx > 2) { cat("\n>>> FAILED -- inspect the 10 worst sides above\n")
  stop("mapping rejected") }
cat("\n>>> bp -> x REPRODUCES GENESPACE's attachments to within", signif(mx, 3), "genes.\n")
cat("written:", OUT_GENE, "|", OUT_CHR, "\n")
