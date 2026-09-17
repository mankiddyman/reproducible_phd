# ---------------------------------------------------------------------------
# DR18a -- probe GENESPACE's (gene order) -> (plot x) mapping.
# Read-only on all inputs. Writes DR/out/DR18a_chr_ord2x.csv only.
# base R + data.table only (genespace module = R 4.2.0, no dplyr)
# ---------------------------------------------------------------------------
suppressMessages(library(data.table))
options(width = 210)

RDA     <- "../genespace/riparian/Nepenthes_gracilis_geneOrder_rSourceData.rda"
COMBBED <- "../genespace/results/combBed.txt"
OUTCSV  <- "DR/out/DR18a_chr_ord2x.csv"
stopifnot(file.exists(RDA), file.exists(COMBBED))

e <- new.env(); load(RDA, envir = e)
P     <- e$srcd$ggplotObj
blk   <- as.data.table(e$srcd$sourceData$blocks)
chrom <- as.data.table(e$srcd$sourceData$chromosomes)
br    <- as.data.table(P$layers[[1]]$data)

# --- 1. STRUCTURES ---------------------------------------------------------
cat("\n================ 1. STRUCTURES ================\n")
cat("\n-- blocks:", nrow(blk), "rows\n");      print(names(blk)); print(head(blk, 3))
cat("\n-- chromosomes:", nrow(chrom), "rows\n"); print(names(chrom)); print(head(chrom, 5))
cat("\n-- braid layer:", nrow(br), "rows\n");   print(names(br));  print(head(br, 3))

# --- 2. CHROMOSOME BOX GEOMETRY --------------------------------------------
cat("\n================ 2. CHROMOSOME BOXES ================\n")
cat("x1 < x2 for every chromosome? ", all(chrom$x1 < chrom$x2), "\n")
chrom[, width := x2 - x1]
if ("length" %in% names(chrom)) {
  cat("cor(width, length) =", cor(chrom$width, chrom$length), "\n")
  cat("width/length ratio range:", paste(signif(range(chrom$width / chrom$length), 8), collapse=" .. "), "\n")
  cat("  (a constant ratio => ONE global scale factor, orientation carried elsewhere)\n")
}
print(chrom[order(genome, chr)][1:12])

# --- 3. BRAID POLYGON ANATOMY ----------------------------------------------
cat("\n================ 3. BRAID ANATOMY ================\n")
stopifnot("blkID" %in% names(br), "x" %in% names(br), "y" %in% names(br))
cat("vertices per blkID:\n"); print(summary(br[, .N, by = blkID]$N))
ylev <- sort(unique(round(c(chrom$y1, chrom$y2), 6)))
cat("chromosome-box y levels:", paste(ylev, collapse = ", "), "\n")
cat("braid y range:", paste(signif(range(br$y), 6), collapse = " .. "), "\n")
eg <- br[blkID == br$blkID[1]][order(y, x)]
cat("\nexample polygon (", eg$blkID[1], ") -- extreme vertices:\n", sep = "")
print(rbind(head(eg, 4), tail(eg, 4)))

# --- 4. ENDPOINT EXTRACTION ------------------------------------------------
ytol <- 1e-6 * diff(range(br$y))
ends <- br[, {
  ymx <- max(y); ymn <- min(y)
  tx <- x[abs(y - ymx) <= ytol]; bx <- x[abs(y - ymn) <= ytol]
  .(top_lo = min(tx), top_hi = max(tx), bot_lo = min(bx), bot_hi = max(bx),
    n_top = length(tx), n_bot = length(bx), y_top = ymx, y_bot = ymn)
}, by = blkID]
cat("\nvertices sitting exactly on the top attachment: "); print(table(ends$n_top))
cat("vertices sitting exactly on the bottom attachment: "); print(table(ends$n_bot))
cat("blocks with an endpoint span of zero width: ",
    sum(ends$top_hi - ends$top_lo == 0 | ends$bot_hi - ends$bot_lo == 0), "\n")

# --- 5. AFFINE FIT, PER CHROMOSOME -----------------------------------------
cat("\n================ 5. AFFINE FIT ================\n")
gcols <- grep("^genome[12]$", names(blk), value = TRUE)
ccols <- grep("^chr[12]$",    names(blk), value = TRUE)
ocols <- grep("^(start|end)Ord[12]$", names(blk), value = TRUE)
idcol <- names(blk)[vapply(blk, function(z) is.character(z) && any(z %in% ends$blkID), logical(1))]
cat("detected -- genome:", paste(gcols, collapse=","), "| chr:", paste(ccols, collapse=","),
    "| ord:", paste(ocols, collapse=","), "| blkID:", paste(idcol, collapse=","), "\n")
if (length(gcols) < 2 || length(ccols) < 2 || length(ocols) < 4 || length(idcol) < 1)
  stop("expected column names not present in blocks -- inspect section 1 printout")
idcol <- idcol[1]

gy  <- chrom[, .(ymid = mean(c(y1, y2))), by = genome]
pts <- rbind(
  blk[, .(blkID = get(idcol), genome = get(gcols[1]), chr = get(ccols[1]),
          ordA = startOrd1, ordB = endOrd1, side = 1L)],
  blk[, .(blkID = get(idcol), genome = get(gcols[2]), chr = get(ccols[2]),
          ordA = startOrd2, ordB = endOrd2, side = 2L)])
pts <- merge(pts, gy,   by = "genome", all.x = TRUE)
pts <- merge(pts, ends, by = "blkID",  all.x = TRUE)
cat("block sides:", nrow(pts), "| with a matched polygon:", sum(!is.na(pts$top_lo)), "\n")
pts <- pts[!is.na(top_lo)]
pts[, is_top := ymid == max(ymid), by = blkID]
pts[, `:=`(xlo = fifelse(is_top, top_lo, bot_lo),
           xhi = fifelse(is_top, top_hi, bot_hi))]

fitchr <- function(d) {
  o  <- c(pmin(d$ordA, d$ordB), pmax(d$ordA, d$ordB))
  if (length(unique(o)) < 2) return(data.table(n=nrow(d), orient=NA_integer_,
      r2=NA_real_, intercept=NA_real_, slope=NA_real_, max_resid=NA_real_))
  go <- function(x) { m <- lm(x ~ o); list(r2=summary(m)$r.squared, a=unname(coef(m)[1]),
                                           b=unname(coef(m)[2]), mx=max(abs(residuals(m)))) }
  fp <- go(c(d$xlo, d$xhi))      # forward: low ord -> low x
  fn <- go(c(d$xhi, d$xlo))      # flipped: low ord -> high x
  bp <- fp$r2 >= fn$r2
  b  <- if (bp) fp else fn
  data.table(n = nrow(d), orient = if (bp) 1L else -1L,
             r2 = b$r2, intercept = b$a, slope = b$b, max_resid = b$mx)
}
res <- pts[, fitchr(.SD), by = .(genome, chr),
           .SDcols = c("ordA","ordB","xlo","xhi")]
res <- merge(res, chrom[, .(genome, chr, box_x1 = x1, box_x2 = x2)],
             by = c("genome","chr"), all.x = TRUE)
print(res[order(r2)], nrows = 100)

cat("\n-- slope magnitudes (one global scale?) --\n"); print(summary(abs(res$slope)))
cat("-- orientation --\n"); print(table(res$orient, useNA = "ifany"))

# predicted x must land inside the chromosome box
pts2 <- merge(pts, res[, .(genome, chr, intercept, slope)], by = c("genome","chr"))
pts2[, `:=`(px1 = intercept + slope * ordA, px2 = intercept + slope * ordB)]
ovr <- merge(pts2[, .(pmin = min(c(px1,px2)), pmax = max(c(px1,px2))), by = .(genome, chr)],
             chrom[, .(genome, chr, x1, x2)], by = c("genome","chr"))
ovr[, overshoot := pmax(x1 - pmin, pmax - x2)]
cat("\nmax overshoot beyond the chromosome box (plot-x units):",
    signif(max(ovr$overshoot, na.rm = TRUE), 5), "\n")
print(ovr[order(-overshoot)][1:5])

# --- 6. combBed: is Mb -> ord exact? ---------------------------------------
cat("\n================ 6. combBed ================\n")
cb <- fread(COMBBED)
print(names(cb)); print(head(cb, 3)); cat("rows:", nrow(cb), "\n")
scol <- grep("^start$", names(cb), value = TRUE)
ocol <- grep("^ord$",   names(cb), value = TRUE)
if (length(scol) && length(ocol)) {
  mono <- cb[, .(rho = if (.N > 3) suppressWarnings(
      cor(get(scol), get(ocol), method = "spearman")) else NA_real_), by = .(genome, chr)]
  cat("Spearman(start, ord) within chromosome, range:",
      paste(signif(range(mono$rho, na.rm = TRUE), 6), collapse = " .. "), "\n")
  cat("chromosomes with rho < 0.999:", sum(mono$rho < 0.999, na.rm = TRUE), "\n")
  if (any(mono$rho < 0.999, na.rm = TRUE)) print(mono[rho < 0.999][order(rho)])
} else cat("!! no start/ord columns -- see names above\n")

# --- 7. VERDICT ------------------------------------------------------------
fwrite(res, OUTCSV)
bad <- sum(res$r2 < 0.9999, na.rm = TRUE)
cat("\n################ VERDICT ################\n")
cat("chromosomes fitted:", nrow(res), "\n")
cat("worst R^2:", format(min(res$r2, na.rm = TRUE), digits = 12), "\n")
cat("chromosomes with R^2 < 0.9999:", bad, "\n")
cat("largest residual, any chromosome (plot-x units):",
    signif(max(res$max_resid, na.rm = TRUE), 5), "\n")
if (is.na(bad) || bad > 0) {
  cat("\n>>> AFFINE HYPOTHESIS REJECTED for", bad, "chromosomes.\n")
  cat(">>> Fall back to polygon warping. Offending chromosomes listed above.\n")
} else {
  cat("\n>>> AFFINE HYPOTHESIS HOLDS: ord -> x is exactly linear per chromosome.\n")
  cat(">>> mb2x can be made exact; no polygon surgery needed.\n")
}
cat("written:", OUTCSV, "\n")
