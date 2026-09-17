#!/usr/bin/env Rscript
# DR21 -- capensis chr15_collapsed in BOTH coordinate systems, side by side.
# Answers one question: is the tract layout displaced, or is the axis rank?
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({ library(data.table) })
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=200)
GS <- file.path(dirname(getwd()), "genespace")
KEY <- "Drosera_capensis chr15_collapsed"

cat("\n== riparian source objects available ==\n")
print(list.files(file.path(GS,"riparian"), pattern="rSourceData", full.names=FALSE))

e <- new.env(); load(file.path(GS,"riparian","Nepenthes_gracilis_geneOrder_rSourceData.rda"), envir=e)
BLK   <- as.data.frame(e$srcd$sourceData$blocks)
CHROM <- as.data.frame(e$srcd$sourceData$chromosomes); CHROM$key <- paste(CHROM$genome, CHROM$chr)
BED <- as.data.frame(fread(file.path(GS,"results","combBed.txt")))
BED$key <- paste(BED$genome, BED$chr); BED <- BED[BED$key %in% CHROM$key, ]
BED <- BED[order(BED$key, BED$ord), ]
BED$r <- ave(BED$ord, BED$key, FUN=function(z) z - min(z) + 1L)
AN <- unique(rbind(
  data.frame(key=paste(BLK$genome1,BLK$chr1), gene=BLK$firstGene1, pord=BLK$startOrd1),
  data.frame(key=paste(BLK$genome1,BLK$chr1), gene=BLK$lastGene1,  pord=BLK$endOrd1),
  data.frame(key=paste(BLK$genome2,BLK$chr2), gene=BLK$firstGene2, pord=BLK$startOrd2),
  data.frame(key=paste(BLK$genome2,BLK$chr2), gene=BLK$lastGene2,  pord=BLK$endOrd2), stringsAsFactors=FALSE))
AN$r <- BED$r[match(paste(AN$key,AN$gene), paste(BED$key,BED$ofID))]; AN <- AN[!is.na(AN$r), ]
OFF <- tapply(AN$r - AN$pord, AN$key, function(z) z[1])
BED$cx1  <- CHROM$x1[match(BED$key,CHROM$key)]; BED$cx2 <- CHROM$x2[match(BED$key,CHROM$key)]
BED$clen <- CHROM$length[match(BED$key,CHROM$key)]
BED$px   <- BED$cx1 + pmin(pmax(BED$r - OFF[BED$key], 1), BED$clen)
BED$mb   <- (BED$start + BED$end)/2e6

G <- read.csv("DR/out/DR02_propagated_blocks.csv", stringsAsFactors=FALSE)
G$key <- paste(G$genome, G$chr); G <- G[G$label %in% c("A","B"), ]
m <- match(paste(BED$key,BED$id), paste(G$key,G$gene))
BED$lab <- G$label[m]; BED$reg <- G$region[m]

z <- BED[BED$key==KEY & !is.na(BED$lab), ]
z <- z[order(z$mb), ]
cat(sprintf("\n== %s ==\n  genes %d | span %.1f-%.1f Mb | box x %.0f..%.0f (%d plot units)\n",
    KEY, nrow(z), min(z$mb), max(z$mb), z$cx1[1], z$cx2[1], z$clen[1]))

cat("\n== gene density by Mb decile (this is the whole story if it is skewed) ==\n")
br <- seq(min(z$mb), max(z$mb), length.out=11)
d  <- data.frame(mb_from=round(head(br,-1),1), mb_to=round(tail(br,-1),1),
                 genes=as.integer(table(cut(z$mb, br, include.lowest=TRUE))))
d$pct_of_genes <- round(100*d$genes/sum(d$genes),1)
print(d)

cat("\n== TRACTS, both coordinate systems ==\n")
r <- rle(paste(z$reg, z$lab)); en <- cumsum(r$lengths); st <- en - r$lengths + 1
TR <- do.call(rbind, lapply(seq_along(r$lengths), function(j){ w <- z[st[j]:en[j], ]
  data.frame(region=w$reg[1], label=w$lab[1], n=nrow(w),
             mb_lo=min(w$mb), mb_hi=max(w$mb), x_lo=min(w$px), x_hi=max(w$px)) }))
TR <- TR[TR$n >= 10, ]
mbspan <- max(z$mb)-min(z$mb); xspan <- z$cx2[1]-z$cx1[1]
TR$mb_pct_start <- round(100*(TR$mb_lo-min(z$mb))/mbspan, 1)
TR$mb_pct_end   <- round(100*(TR$mb_hi-min(z$mb))/mbspan, 1)
TR$x_pct_start  <- round(100*(TR$x_lo-z$cx1[1])/xspan, 1)
TR$x_pct_end    <- round(100*(TR$x_hi-z$cx1[1])/xspan, 1)
TR$genes_per_mb <- round(TR$n/pmax(TR$mb_hi-TR$mb_lo, 0.01))
print(TR[, c("region","label","n","mb_lo","mb_hi","mb_pct_start","mb_pct_end",
             "x_pct_start","x_pct_end","genes_per_mb")], row.names=FALSE)
cat(sprintf("\n  tract order identical in both systems? %s\n",
    identical(order(TR$mb_lo), order(TR$x_lo))))
cat("  IF ORDER IS IDENTICAL BUT THE PERCENTAGES DIFFER, THE FIGURE IS CORRECT\n")
cat("  AND THE DISPLACEMENT YOU SEE IS THE RANK AXIS, NOT A BUG.\n")

cat("\n== BRAID ATTACHMENTS on this chromosome, both systems ==\n")
S <- rbind(
  data.frame(blkID=BLK$blkID, key=paste(BLK$genome1,BLK$chr1), other=paste(BLK$genome2,BLK$chr2),
             bp1=BLK$startBp1, bp2=BLK$endBp1, g1=BLK$firstGene1, g2=BLK$lastGene1, stringsAsFactors=FALSE),
  data.frame(blkID=BLK$blkID, key=paste(BLK$genome2,BLK$chr2), other=paste(BLK$genome1,BLK$chr1),
             bp1=BLK$startBp2, bp2=BLK$endBp2, g1=BLK$firstGene2, g2=BLK$lastGene2, stringsAsFactors=FALSE))
S <- S[S$key==KEY, ]
kk <- paste(BED$key, BED$ofID)
S$xa <- BED$px[match(paste(S$key,S$g1), kk)]; S$xb <- BED$px[match(paste(S$key,S$g2), kk)]
S$mba <- pmin(S$bp1,S$bp2)/1e6; S$mbb <- pmax(S$bp1,S$bp2)/1e6
S$mb_pct  <- round(100*(S$mba-min(z$mb))/mbspan, 1)
S$x_pct   <- round(100*(pmin(S$xa,S$xb)-z$cx1[1])/xspan, 1)
tract_at <- function(v, lo, hi) { k <- which(lo <= v & hi >= v); if (length(k)) paste0(TR$region[k[1]],"/",TR$label[k[1]]) else "-" }
S$tract_mb <- mapply(tract_at, (S$mba+S$mbb)/2, MoreArgs=list(lo=TR$mb_lo, hi=TR$mb_hi))
S$tract_x  <- mapply(tract_at, (S$xa+S$xb)/2,  MoreArgs=list(lo=TR$x_lo,  hi=TR$x_hi))
S$agree <- S$tract_mb == S$tract_x
print(S[order(S$mba), c("other","mba","mbb","mb_pct","x_pct","tract_mb","tract_x","agree")], row.names=FALSE)
cat(sprintf("\n  attachments landing in the SAME tract in both systems: %d of %d\n",
            sum(S$agree), nrow(S)))
if (all(S$agree)) cat("  => the data agree. The figure differs from the zoom ONLY by axis scaling.\n") else
  cat("  => REAL BUG: some attachments change tract between bp and rank. Rows above with agree=FALSE.\n")
