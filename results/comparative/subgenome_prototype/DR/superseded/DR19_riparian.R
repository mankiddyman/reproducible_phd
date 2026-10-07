#!/usr/bin/env Rscript
# ============================================================================
# DR19 — DR17, with ONE canonical gene -> plot x map
#
# DR18 established: GENESPACE plots a SUBSET of combBed genes, and
#     plot_ord = (rank within chromosome by ord) - offset(chromosome)
#     x        = chromosome$x1 + plot_ord
# offset is a constant per chromosome (0..98), verified across 4184 anchor
# genes with zero spread. No chromosome is flipped.
#
# DR17 computed gene->x in three different places, all wrong:
#   bars      : rank without the offset
#   endlab    : rank rescaled by combBed count over box width
#   Dionaea   : same rescale
# Here there is exactly one map, BED$px, used everywhere, and it is proved
# against GENESPACE's own polygon attachments by GENE IDENTITY before drawing.
#
# IN   ../genespace/riparian/Nepenthes_gracilis_geneOrder_rSourceData.rda
#      ../genespace/results/combBed.txt
#      DR/out/DR02_propagated_blocks.csv, fractionation_by_chrpair.csv
# OUT  DR/fig/DR19_riparian_AB.{pdf,png}, DR/out/DR19_bar_braid_mismatch.csv
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({ library(GENESPACE); library(data.table); library(ggplot2) })
setwd(Sys.getenv("SUBG_BASE", getwd()))
GS <- file.path(dirname(getwd()), "genespace")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }

hr("1. LOAD THE SAVED RIPARIAN")
e <- new.env()
load(file.path(GS,"riparian","Nepenthes_gracilis_geneOrder_rSourceData.rda"), envir=e)
P     <- e$srcd$ggplotObj
BR    <- as.data.frame(P$layers[[1]]$data)
CHROM <- as.data.frame(e$srcd$sourceData$chromosomes)
BLK   <- as.data.frame(e$srcd$sourceData$blocks)
CHROM$key <- paste(CHROM$genome, CHROM$chr)
cat(sprintf("  braid rows %d | blocks %d | chromosomes %d\n",
            nrow(BR), nrow(BLK), nrow(CHROM)))

hr("2. CANONICAL GENE -> PLOT X")
BED <- as.data.frame(fread(file.path(GS,"results","combBed.txt")))
BED$key <- paste(BED$genome, BED$chr)
BED <- BED[BED$key %in% CHROM$key, ]
BED <- BED[order(BED$key, BED$ord), ]
BED$r <- ave(BED$ord, BED$key, FUN = function(z) z - min(z) + 1L)
cat(sprintf("  combBed genes on plotted chromosomes: %d\n", nrow(BED)))

AN <- unique(rbind(
  data.frame(key=paste(BLK$genome1,BLK$chr1), gene=BLK$firstGene1, pord=BLK$startOrd1),
  data.frame(key=paste(BLK$genome1,BLK$chr1), gene=BLK$lastGene1,  pord=BLK$endOrd1),
  data.frame(key=paste(BLK$genome2,BLK$chr2), gene=BLK$firstGene2, pord=BLK$startOrd2),
  data.frame(key=paste(BLK$genome2,BLK$chr2), gene=BLK$lastGene2,  pord=BLK$endOrd2),
  stringsAsFactors=FALSE))
AN$r <- BED$r[match(paste(AN$key, AN$gene), paste(BED$key, BED$ofID))]
AN <- AN[!is.na(AN$r), ]
AN$d <- AN$r - AN$pord
spread <- tapply(AN$d, AN$key, function(z) diff(range(z)))
cat(sprintf("  anchors %d on %d chromosomes | offset spread max %g\n",
            nrow(AN), length(spread), max(spread)))
if (max(spread) != 0) die("offset is not constant on %d chromosomes",
                          sum(spread != 0))
OFF <- tapply(AN$d, AN$key, function(z) z[1])
if (!all(CHROM$key %in% names(OFF))) die("no anchors on %d chromosomes",
                                         sum(!CHROM$key %in% names(OFF)))
cat(sprintf("  offsets: min %g median %g max %g\n",
            min(OFF), median(OFF), max(OFF)))

BED$off      <- OFF[BED$key]
BED$cx1      <- CHROM$x1[match(BED$key, CHROM$key)]
BED$clen     <- CHROM$length[match(BED$key, CHROM$key)]
BED$plot_ord <- BED$r - BED$off
BED$px       <- BED$cx1 + pmin(pmax(BED$plot_ord, 1), BED$clen)
cat(sprintf("  genes outside the plotted range, clamped to the box edge: %d\n",
            sum(BED$plot_ord < 1 | BED$plot_ord > BED$clen)))

hr("2b. PROOF — predicted x vs GENESPACE's own polygon attachments")
EP <- as.data.frame(as.data.table(BR)[, .(xlo=min(x), xhi=max(x)),
                                      by=.(blkID, yrow=round(y,6))])
SD <- rbind(
  data.frame(blkID=BLK$blkID, yrow=round(as.numeric(BLK$y1),6),
             key=paste(BLK$genome1,BLK$chr1), ga=BLK$firstGene1, gb=BLK$lastGene1,
             stringsAsFactors=FALSE),
  data.frame(blkID=BLK$blkID, yrow=round(as.numeric(BLK$y2),6),
             key=paste(BLK$genome2,BLK$chr2), ga=BLK$firstGene2, gb=BLK$lastGene2,
             stringsAsFactors=FALSE))
SD <- merge(SD, EP, by=c("blkID","yrow"))
kk <- paste(BED$key, BED$ofID)
SD$xa <- BED$px[match(paste(SD$key, SD$ga), kk)]
SD$xb <- BED$px[match(paste(SD$key, SD$gb), kk)]
SD <- SD[!is.na(SD$xa) & !is.na(SD$xb), ]
SD$err <- pmax(abs(pmin(SD$xa,SD$xb) - SD$xlo), abs(pmax(SD$xa,SD$xb) - SD$xhi))
cat(sprintf("  block sides checked: %d | max error: %.3g plot units\n",
            nrow(SD), max(SD$err)))
if (nrow(SD) < 2000) die("only %d block sides resolvable", nrow(SD))
if (max(SD$err) > 1e-6) {
  print(head(SD[order(-SD$err), c("key","ga","gb","xa","xb","xlo","xhi","err")], 10))
  die("gene -> x does not reproduce GENESPACE's attachments")
}
cat("  EXACT. The bars and the braids are now in one coordinate system.\n")

hr("3. GENES, LABELS, TRACTS")
BR$base <- sub("_[0-9]+ chr[^ ]*$", "", BR$blkID)
G <- read.csv("DR/out/DR02_propagated_blocks.csv", stringsAsFactors=FALSE)
cat("  propagated genes:", nrow(G), "| columns:", paste(names(G), collapse=", "), "\n")
gcol <- intersect(c("gene","id","ofID"), names(G))[1]
posc <- intersect(c("mid","start","pos","bp"), names(G))[1]
G$mb <- G[[posc]] / 1e6
for (nm in c("label","region")) if (!nm %in% names(G)) {
  SS <- read.csv("DR/out/DR02_segments.csv", stringsAsFactors=FALSE)
  G[[nm]] <- SS[[nm]][match(G$segment, SS$segment)] }
G$lab <- G$label
G <- G[!is.na(G$lab) & G$lab %in% c("A","B") & !is.na(G$mb), ]

# exact join first: DR02 gene id -> combBed
G$key <- paste(G$genome, G$chr)
m_id  <- match(paste(G$key, G[[gcol]]), paste(BED$key, BED$id))
m_of  <- match(paste(G$key, G[[gcol]]), kk)
use   <- if (sum(!is.na(m_id)) >= sum(!is.na(m_of))) m_id else m_of
cat(sprintf("  exact id join: %d of %d genes (%.1f%%)\n",
            sum(!is.na(use)), nrow(G), 100*mean(!is.na(use))))
G$px <- BED$px[use]
# fallback for stragglers: interpolate on gene midpoint within the chromosome
nmiss <- sum(is.na(G$px))
if (nmiss) {
  BED$bmb <- (BED$start + BED$end)/2e6
  for (k in unique(G$key[is.na(G$px)])) {
    b <- BED[BED$key == k, ]; if (!nrow(b)) next
    b <- b[order(b$bmb), ]
    i <- which(G$key == k & is.na(G$px))
    G$px[i] <- approx(b$bmb, b$px, xout=G$mb[i], rule=2)$y } }
cat(sprintf("  interpolated fallback: %d | unplaced (dropped): %d\n",
            nmiss - sum(is.na(G$px)), sum(is.na(G$px))))
G <- G[!is.na(G$px), ]

CHR <- data.frame(genome=CHROM$genome, chr=CHROM$chr, x0=CHROM$x1, x1=CHROM$x2,
                  y=(CHROM$y1 + CHROM$y2)/2, key=CHROM$key, stringsAsFactors=FALSE)
G$yrow <- CHR$y[match(G$key, CHR$key)]
G <- G[!is.na(G$yrow), ]
cat(sprintf("  genes placed in plot coordinates: %d\n", nrow(G)))

# tracts: contiguous runs of (region,label) along the chromosome, in PLOT x
TR <- do.call(rbind, lapply(split(G, G$key), function(z){
  z <- z[order(z$px), ]
  r <- rle(paste(z$region, z$lab)); en <- cumsum(r$lengths); st <- en - r$lengths + 1
  do.call(rbind, lapply(seq_along(r$lengths), function(j){
    w <- z[st[j]:en[j], ]
    data.frame(genome=w$genome[1], region=w$region[1], chr=w$chr[1], label=w$lab[1],
               xa=min(w$px), xb=max(w$px), n=nrow(w), stringsAsFactors=FALSE) }))
}))
TR <- TR[TR$n >= 10, ]
.r5 <- TR[TR$genome=="Drosera_regia" & TR$chr=="chr5_collapsed", ]
if (length(unique(.r5$label)) != 2) { print(.r5); die("regia chr5_collapsed lost its A/B split") }
.both <- sum(vapply(split(TR, paste(TR$genome,TR$chr)),
                    function(z) length(unique(z$label)) > 1, logical(1)))
cat(sprintf("  tracts: %d (%d A, %d B) | chromosomes with BOTH labels: %d\n",
            nrow(TR), sum(TR$label=="A"), sum(TR$label=="B"), .both))
if (.both < 25) die("only %d chromosomes carry both labels -- tracts merged", .both)
TR$y <- CHR$y[match(paste(TR$genome,TR$chr), paste(CHR$genome,CHR$chr))]
TR <- TR[!is.na(TR$y), ]

hr("4. DIONAEA + INTERVALS")
PAIRS <- read.csv("fractionation_by_chrpair.csv", stringsAsFactors=FALSE)
s1 <- PAIRS$retained_more == PAIRS$chrA
DIOSIDE <- c(setNames(ifelse(s1,"A","B"), PAIRS$chrA),
             setNames(ifelse(s1,"B","A"), PAIRS$chrB))
DC <- CHR[CHR$genome=="Dionaea_muscipula", ]
DC$label <- unname(DIOSIDE[DC$chr])
if (any(is.na(DC$label))) die("Dionaea chromosomes without an A/B side: %s",
                              paste(DC$chr[is.na(DC$label)], collapse=", "))
IV <- rbind(TR[, c("genome","chr","region","label","xa","xb","y")],
            data.frame(genome=DC$genome, chr=DC$chr, region="Dionaea",
                       label=DC$label, xa=DC$x0, xb=DC$x1, y=DC$y, stringsAsFactors=FALSE))
IV$ytop <- IV$y - 0.05; IV$ybot <- IV$y - 0.11
cat(sprintf("  intervals to draw: %d\n", nrow(IV)))

DIOG <- BED[BED$genome=="Dionaea_muscipula", ]
DIOG$lab <- unname(DIOSIDE[DIOG$chr])
DIOG <- DIOG[!is.na(DIOG$lab), ]
DIOG$yrow <- CHR$y[match(DIOG$key, CHR$key)]
ALLG <- rbind(data.frame(x=G$px, yrow=G$yrow, lab=G$lab, reg=G$region, stringsAsFactors=FALSE),
              data.frame(x=DIOG$px, yrow=DIOG$yrow, lab=DIOG$lab, reg=NA_character_,
                         stringsAsFactors=FALSE))
ALLG <- ALLG[!is.na(ALLG$x) & !is.na(ALLG$yrow), ]
cat(sprintf("  labelled genes on the plot: %d\n", nrow(ALLG)))

hr("5. LABEL THE BRAID ENDS")
endlab <- function(x0, x1, y, reg) {
  k <- which(abs(ALLG$yrow - y) < 0.25 & ALLG$x >= x0 & ALLG$x <= x1 &
             (is.na(ALLG$reg) | ALLG$reg == reg))
  if (length(k) < 3) return(NA_character_)
  t <- table(ALLG$lab[k]); if (max(t)/sum(t) < 0.95) return(NA_character_)
  names(t)[which.max(t)] }
BRs <- split(seq_len(nrow(BR)), BR$blkID)
res <- do.call(rbind, lapply(names(BRs), function(id) {
  ix <- BRs[[id]]; yy <- BR$y[ix]; xx <- BR$x[ix]
  lo <- which(yy < min(yy)+0.02); hi <- which(yy > max(yy)-0.02)
  if (!length(lo) || !length(hi)) return(NULL)
  rg <- sub("^.*_[0-9]+ ", "", id)
  data.frame(blkID=id, region=rg, y1=min(yy), y2=max(yy),
             x1a=min(xx[lo]), x1b=max(xx[lo]), x2a=min(xx[hi]), x2b=max(xx[hi]),
             L1=endlab(min(xx[lo]), max(xx[lo]), min(yy), rg),
             L2=endlab(min(xx[hi]), max(xx[hi]), max(yy), rg), stringsAsFactors=FALSE) }))
nep_y <- CHR$y[CHR$genome=="Nepenthes_gracilis"][1]
res$nep1 <- abs(res$y1-nep_y) < 0.25; res$nep2 <- abs(res$y2-nep_y) < 0.25
res$cls <- ifelse(!is.na(res$L1) & !is.na(res$L2),
             ifelse(res$L1==res$L2, res$L1, "discordant"),
           ifelse(is.na(res$L1) & res$nep1 & !is.na(res$L2), res$L2,
           ifelse(is.na(res$L2) & res$nep2 & !is.na(res$L1), res$L1, "no call")))
cat("  braids by class:\n"); print(table(res$cls))
BR$cls <- res$cls[match(BR$blkID, res$blkID)]; BR$cls[is.na(BR$cls)] <- "no call"

hr("6. AUDIT — every drawn braid end against the bar beneath it")
barlab_at <- function(x, y, reg) {
  k <- which(abs(IV$y - y) < 0.25 & IV$xa <= x & IV$xb >= x)
  if (!length(k)) return(NA_character_)
  ks <- k[IV$region[k] == reg | IV$region[k] == "Dionaea"]
  if (length(ks)) IV$label[ks[1]] else NA_character_ }
D <- res[res$cls %in% c("A","B"), ]
D$bar1 <- mapply(barlab_at, (D$x1a+D$x1b)/2, D$y1, D$region)
D$bar2 <- mapply(barlab_at, (D$x2a+D$x2b)/2, D$y2, D$region)
D$mm1 <- !is.na(D$bar1) & D$bar1 != D$cls
D$mm2 <- !is.na(D$bar2) & D$bar2 != D$cls
MM <- D[D$mm1 | D$mm2, ]
cat(sprintf("  drawn braids: %d\n", nrow(D)))
cat(sprintf("  end 1 on a wrong-coloured bar: %d (%.2f%%)\n", sum(D$mm1), 100*mean(D$mm1)))
cat(sprintf("  end 2 on a wrong-coloured bar: %d (%.2f%%)\n", sum(D$mm2), 100*mean(D$mm2)))
if (nrow(MM)) {
  cat("\n  by region:\n"); print(head(as.data.frame(sort(table(MM$region), decreasing=TRUE)), 8))
  write.csv(MM[, c("blkID","region","cls","bar1","bar2","y1","y2")],
            "DR/out/DR19_bar_braid_mismatch.csv", row.names=FALSE)
  cat("  wrote DR/out/DR19_bar_braid_mismatch.csv\n")
} else cat("\n  ZERO mismatches -- bars and braids agree everywhere\n")

hr("7. DRAW")
lum <- function(h_,s_,v_){h <- rgb2hsv(col2rgb(h_)); hsv(h["h",],pmin(1,h["s",]*s_),pmin(1,h["v",]*v_))}
cols <- unique(BR$color)
Ac <- setNames(lum(cols,0.55,1.10), cols); Bc <- setNames(lum(cols,1.00,0.88), cols)
BR$fillcol <- ifelse(BR$cls=="A", Ac[BR$color], Bc[BR$color])
IV$barcol  <- ifelse(IV$label=="A", "#1D9E75", "#D85A30")
PAL <- c(setNames(unique(BR$fillcol),unique(BR$fillcol)),
         setNames(unique(IV$barcol), unique(IV$barcol)))
mkl <- function(d,al) layer(geom="polygon", stat="identity", position="identity",
  data=d, mapping=aes(x=x,y=y,group=blkID,fill=fillcol),
  params=list(alpha=al,colour=NA,linewidth=0), show.legend=FALSE)
p <- P
p$layers <- c(list(mkl(BR[BR$cls=="A",],0.55), mkl(BR[BR$cls=="B",],0.88)),
              P$layers[-1],
              list(layer(geom="rect", stat="identity", position="identity", data=IV,
                mapping=aes(xmin=xa,xmax=xb,ymin=ybot,ymax=ytop,fill=barcol),
                params=list(colour=NA), show.legend=FALSE)))
p$scales$scales <- p$scales$scales[
  !vapply(p$scales$scales, function(sc) "fill" %in% sc$aesthetics, logical(1))]
p <- p + scale_fill_manual(values=PAL, guide="none") + theme_minimal(11) +
  theme(panel.background=element_rect(fill="grey96",colour=NA),
        plot.background=element_rect(fill="grey96",colour=NA),
        panel.grid=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(),
        axis.text.y=element_text(face="italic",size=11),
        plot.title=element_text(face="bold",size=14),
        plot.subtitle=element_text(size=9,colour="grey30")) +
  labs(title="Droseraceae riparian, phased by subgenome",
       subtitle=paste0("Bars beneath each chromosome are the phased tracts: ",
         "green = A, orange = B. Braids are filtered against the same\n",
         "tracts -- drawn only where both ends fall on one subgenome. ",
         "MUTED = A, VIVID = B. Hue = ancestral region."),
       x="Chromosomes scaled by gene rank order", y=NULL)
ggsave("DR/fig/DR19_riparian_AB.pdf", p, width=16, height=9)
ggsave("DR/fig/DR19_riparian_AB.png", p, width=16, height=9, dpi=200)
cat("  wrote DR/fig/DR19_riparian_AB.{pdf,png}\n")
cat(sprintf("\n  FINAL: %d drawn braids, %d with an end on a wrong-coloured bar\n",
            nrow(D), nrow(MM)))
