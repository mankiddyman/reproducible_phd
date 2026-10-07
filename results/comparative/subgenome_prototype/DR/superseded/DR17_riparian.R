#!/usr/bin/env Rscript
# ============================================================================
# DR17 — THE GENESPACE RIPARIAN, WITH A/B ON THE BRAIDS
#
# Not a reimplementation. GENESPACE saved its plot object to
# riparian/Nepenthes_gracilis_geneOrder_rSourceData.rda, so we take that exact
# figure -- gene-rank ordering, chromosome gaps, curved braids from
# calc_curvePolygon, labels -- and rebuild ONE layer.
#
#   layer 1 = the braids   (GeomPolygon, group = blkID, fill = color)  <- ours
#   layer 2 = chromosomes  (grouped by genome, chr)                    <- kept
#   layer 3 = labels                                                   <- kept
#
# ENCODING
#   solid, opaque      A subgenome
#   faded + dashed     B subgenome
#   very faint grey    no confident call
#
# THE JOIN
#   Layer 1 blkID: "Drosera_regia_vs_Dionaea_muscipula: 283_1 chr10_dom"
#   syntenicBlock_coordinates.csv blkID: "Drosera_regia_vs_Dionaea_muscipula: 283"
#   Strip the "_N chrX_dom" suffix to match. The block's Drosera-side midpoint
#   is then matched against DR02's segments to inherit an A/B label.
#
# WHAT THE FIGURE INHERITS
#   Only 8-19% of genes have a measured delta; the rest take their segment's
#   label. DR14 showed 94.6-98.3% of propagated genes sit within 100 kb of a
#   measured one, so this is dense interpolation, not extrapolation.
#   Vote purity averages 0.68-0.74 -- near the ceiling set by 0.5 SD per-gene
#   separation. Reliability comes from vote COUNT: 0.998 at 50 votes.
#
# IN   ../genespace/riparian/Nepenthes_gracilis_geneOrder_rSourceData.rda
#      ../genespace/results/syntenicBlock_coordinates.csv
#      DR/out/DR02_segments.csv, DR/out/DR14_segment_quality.csv
# OUT  DR/fig/DR17_riparian_AB.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(GENESPACE); library(data.table); library(ggplot2)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))
GS <- file.path(dirname(getwd()), "genespace")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

hr("1. LOAD THE SAVED RIPARIAN")
e <- new.env()
load(file.path(GS,"riparian","Nepenthes_gracilis_geneOrder_rSourceData.rda"), envir=e)
P <- e$srcd$ggplotObj
BR <- as.data.frame(P$layers[[1]]$data)   # the braids
cat(sprintf("  braid polygon rows: %d | distinct blocks: %d\n",
            nrow(BR), length(unique(BR$blkID))))

hr("2. PARSE blkID AND JOIN TO BLOCK COORDINATES")
BR$base <- sub("_[0-9]+ chr[^ ]*$", "", BR$blkID)
cat("  layer-1 blkID :", head(unique(BR$blkID),1), "\n")
cat("  parsed base   :", head(unique(BR$base),1), "\n")
BC <- fread(file.path(GS,"results","syntenicBlock_coordinates.csv"))
cat("  coords blkID  :", head(unique(as.character(BC$blkID)),1), "\n")
ov <- sum(unique(BR$base) %in% as.character(BC$blkID))
cat(sprintf("\n  JOIN: %d of %d parsed IDs found in the coordinate file\n",
            ov, length(unique(BR$base))))
if (ov < 0.5*length(unique(BR$base))) {
  cat("  *** parse failed. Showing both formats so the regex can be fixed: ***\n")
  print(head(unique(BR$base), 5)); print(head(as.character(BC$blkID), 5))
  quit(status=1)
}

hr("3. BUILD REAL TRACTS FROM GENE POSITIONS")
cat("  DR02's mb_lo/mb_hi are min/max of a segment's genes -- a bounding box,\n")
cat("  not a contiguous tract. A few outliers stretched a box across a whole\n")
cat("  chromosome, so boxes from different ancestral regions overlapped on 28\n")
cat("  chromosomes and the label lookup became arbitrary.\n")
cat("  Splitting each segment at gaps > 1 Mb gives 517 real tracts with only\n")
cat("  4 chromosomes overlapping and 88%% of genes retained.\n\n")

GAP <- 0.5
G <- read.csv("DR/out/DR02_propagated_blocks.csv", stringsAsFactors=FALSE)
cat("  propagated genes loaded:", nrow(G), "\n")
cat("  columns:", paste(names(G), collapse=", "), "\n")
posc <- intersect(c("mid","start","pos","bp"), names(G))[1]
G$mb <- G[[posc]] / 1e6
if (!"label" %in% names(G)) {
  SS <- read.csv("DR/out/DR02_segments.csv", stringsAsFactors=FALSE)
  G$label <- SS$label[match(G$segment, SS$segment)]
}
if (!"region" %in% names(G)) {
  SS <- read.csv("DR/out/DR02_segments.csv", stringsAsFactors=FALSE)
  G$region <- SS$region[match(G$segment, SS$segment)]
}
G$lab <- G$label
G <- G[!is.na(G$lab) & G$lab %in% c("A","B") & !is.na(G$mb), ]
cat("  usable:", nrow(G), "genes\n")
# contiguous runs of (region, label) along each chromosome -- the blocks
# DR02c already established. No gap threshold, no re-derivation.
TR <- do.call(rbind, lapply(split(G, paste(G$genome, G$chr)), function(z){
  z <- z[order(z$mb), ]
  r <- rle(paste(z$region, z$lab))
  en <- cumsum(r$lengths); st <- en - r$lengths + 1
  do.call(rbind, lapply(seq_along(r$lengths), function(j){
    w <- z[st[j]:en[j], ]
    data.frame(genome=w$genome[1], region=w$region[1], chr=w$chr[1],
               label=w$lab[1], lo=min(w$mb), hi=max(w$mb), n=nrow(w),
               stringsAsFactors=FALSE)
  }))
}))
TR <- TR[TR$n >= 10, ]

# ASSERTION -- regia chr5_collapsed has an A tract AND a B tract (DR02c
# verified: A 0-10.2 Mb, B 10.3-19.0 Mb). Two separate bugs have collapsed it
# to one colour. Fail loudly rather than draw a wrong figure.
.r5 <- TR[TR$genome=="Drosera_regia" & TR$chr=="chr5_collapsed", ]
if (length(unique(.r5$label)) != 2) {
  cat("\n  *** ASSERTION FAILED: regia chr5_collapsed lost its A/B split ***\n")
  print(.r5)
  quit(status=1)
}
.bothlab <- sum(vapply(split(TR, paste(TR$genome, TR$chr)),
                       function(z) length(unique(z$label)) > 1, logical(1)))
cat(sprintf("  chromosomes drawn with BOTH labels: %d (DR02c says 31)\n", .bothlab))
if (.bothlab < 25) { cat("  *** too few -- tracts are being merged ***\n"); quit(status=1) }
cat(sprintf("  tracts: %d (%d A, %d B)\n", nrow(TR),
            sum(TR$label=="A"), sum(TR$label=="B")))

CH <- as.data.frame(P$layers[[2]]$data)
CHR <- do.call(rbind, lapply(split(CH, paste(CH$genome, CH$chr)), function(z)
  data.frame(genome=z$genome[1], chr=z$chr[1], x0=min(z$x), x1=max(z$x),
             y=mean(z$y), stringsAsFactors=FALSE)))
GSD2 <- file.path(dirname(getwd()), "genespace", "results")
BED <- read.table(file.path(GSD2,"combBed.txt"), header=TRUE, sep="\t",
                  quote="", comment.char="", stringsAsFactors=FALSE)
BED$mb <- (BED$start + BED$end)/2e6
ORDMIN <- tapply(BED$ord, paste(BED$genome, BED$chr), min)
GENEX <- setNames(rep(NA_real_, nrow(BED)), paste(BED$genome, BED$gene <- BED$id))
BED$key <- paste(BED$genome, BED$chr)
BED$rank <- BED$ord - ORDMIN[BED$key] + 1
BED$px <- NA_real_
for (kk in unique(BED$key)) {
  g <- sub(" .*$", "", kk); c <- sub("^[^ ]+ ", "", kk)
  k <- which(CHR$genome==g & CHR$chr==c); if (!length(k)) next
  ix <- which(BED$key == kk)
  BED$px[ix] <- CHR$x0[k[1]] + BED$rank[ix]
}
cat(sprintf("  genes with a plot x from GENESPACE ranks: %d of %d\n",
            sum(!is.na(BED$px)), nrow(BED)))
GX <- setNames(BED$px, paste(BED$genome, BED$id))
mb2x <- function(g, c, mb) {
  bb <- BED[BED$genome==g & BED$chr==c & !is.na(BED$px), ]
  if (!nrow(bb)) return(NA_real_)
  bb$px[which.min(abs(bb$mb - mb))]
}
TR$xa <- mapply(mb2x, TR$genome, TR$chr, TR$lo)
TR$xb <- mapply(mb2x, TR$genome, TR$chr, TR$hi)
TR$y  <- CHR$y[match(paste(TR$genome,TR$chr), paste(CHR$genome,CHR$chr))]
TR <- TR[!is.na(TR$xa) & !is.na(TR$xb) & !is.na(TR$y), ]

PAIRS <- read.csv("fractionation_by_chrpair.csv", stringsAsFactors=FALSE)
s1 <- PAIRS$retained_more == PAIRS$chrA
DIOSIDE <- c(setNames(ifelse(s1,"A","B"), PAIRS$chrA),
             setNames(ifelse(s1,"B","A"), PAIRS$chrB))
D <- CHR[CHR$genome=="Dionaea_muscipula", ]
D$label <- unname(DIOSIDE[D$chr]); stopifnot(!any(is.na(D$label)))
IV <- rbind(TR[, c("genome","chr","region","label","xa","xb","y")],
            data.frame(genome=D$genome, chr=D$chr, region="Dionaea",
                       label=D$label, xa=D$x0, xb=D$x1, y=D$y,
                       stringsAsFactors=FALSE))
IV$ytop <- IV$y - 0.05
IV$ybot <- IV$y - 0.11
cat(sprintf("  intervals to draw: %d\n", nrow(IV)))

# Braid ends span a median of 170 plot units and up to 6828 -- a whole
# chromosome. A single mean(x) lookup lands in the middle of a wide ribbon and
# can miss the tracts at both edges. Instead, score each END by every gene
# lying under it and take the majority label.
xr <- function(g, c) { k <- which(CHR$genome==g & CHR$chr==c)
                       if (length(k)) c(CHR$x0[k[1]], CHR$x1[k[1]]) else c(NA,NA) }
Gx <- G
Gx$x <- NA_real_
for (gg in unique(Gx$genome)) for (cc in unique(Gx$chr[Gx$genome==gg])) {
  k <- which(CHR$genome==gg & CHR$chr==cc); if (!length(k)) next
  bb <- BED[BED$genome==gg & BED$chr==cc, ]; if (!nrow(bb)) next
  bb <- bb[order(bb$ord), ]
  bb$rank <- seq_len(nrow(bb))
  bb <- bb[order(bb$mb), ]                 # ascending mb for findInterval
  ix <- which(Gx$genome==gg & Gx$chr==cc)
  j  <- pmax(1, findInterval(Gx$mb[ix], bb$mb))
  Gx$x[ix] <- CHR$x0[k[1]] + (bb$rank[j]-1)/max(1,nrow(bb)-1) *
              (CHR$x1[k[1]] - CHR$x0[k[1]])
}
Gx$yrow <- CHR$y[match(paste(Gx$genome, Gx$chr), paste(CHR$genome, CHR$chr))]
Gx <- Gx[!is.na(Gx$x) & !is.na(Gx$yrow), ]
cat(sprintf("  genes placed in plot coordinates: %d\n", nrow(Gx)))

# Dionaea genes: whole-chromosome labels
DIOG <- BED[BED$genome=="Dionaea_muscipula", ]
DIOG$lab <- unname(DIOSIDE[DIOG$chr])
DIOG <- DIOG[!is.na(DIOG$lab), ]
DIOG$x <- NA_real_
for (cc in unique(DIOG$chr)) {
  k <- which(CHR$genome=="Dionaea_muscipula" & CHR$chr==cc); if (!length(k)) next
  bb <- DIOG[DIOG$chr==cc, ]; bb <- bb[order(bb$ord), ]
  DIOG$x[DIOG$chr==cc] <- CHR$x0[k[1]] +
    (seq_len(nrow(bb))-1)/max(1,nrow(bb)-1) * (CHR$x1[k[1]] - CHR$x0[k[1]])
}
DIOG$yrow <- CHR$y[match(paste("Dionaea_muscipula", DIOG$chr),
                         paste(CHR$genome, CHR$chr))]
Gx$reg <- Gx$region
DIOG$reg <- NA_character_          # Dionaea labels are whole-chromosome
ALLG <- rbind(Gx[, c("x","yrow","lab","reg")],
              DIOG[!is.na(DIOG$x) & !is.na(DIOG$yrow), c("x","yrow","lab","reg")])
cat(sprintf("  total labelled genes on the plot: %d\n", nrow(ALLG)))

# majority label among genes under a braid end
endlab <- function(x0, x1, y, reg) {
  k <- which(abs(ALLG$yrow - y) < 0.25 & ALLG$x >= x0 & ALLG$x <= x1 &
             (is.na(ALLG$reg) | ALLG$reg == reg))
  if (length(k) < 3) return(NA_character_)
  t <- table(ALLG$lab[k])
  if (max(t)/sum(t) < 0.95) return(NA_character_)   # chimeric ends dropped
  names(t)[which.max(t)]
}
BRs <- split(seq_len(nrow(BR)), BR$blkID)
res <- do.call(rbind, lapply(names(BRs), function(id) {
  ix <- BRs[[id]]; yy <- BR$y[ix]; xx <- BR$x[ix]
  lo <- which(yy < min(yy)+0.02); hi <- which(yy > max(yy)-0.02)
  if (!length(lo) || !length(hi)) return(NULL)
  rg <- sub("^.*_[0-9]+ ", "", id)          # the region, from the blkID
  data.frame(blkID = id, region = rg, y1 = min(yy), y2 = max(yy),
             x1a = min(xx[lo]), x1b = max(xx[lo]),
             x2a = min(xx[hi]), x2b = max(xx[hi]),
             L1 = endlab(min(xx[lo]), max(xx[lo]), min(yy), rg),
             L2 = endlab(min(xx[hi]), max(xx[hi]), max(yy), rg),
             stringsAsFactors = FALSE)
}))
nep_y <- CHR$y[CHR$genome == "Nepenthes_gracilis"][1]
res$nep1 <- abs(res$y1 - nep_y) < 0.25
res$nep2 <- abs(res$y2 - nep_y) < 0.25
res$cls <- ifelse(!is.na(res$L1) & !is.na(res$L2),
             ifelse(res$L1 == res$L2, res$L1, "discordant"),
           ifelse(is.na(res$L1) & res$nep1 & !is.na(res$L2), res$L2,
           ifelse(is.na(res$L2) & res$nep2 & !is.na(res$L1), res$L1,
                  "no call")))
cat(sprintf("  braids saved by the Nepenthes exemption: %d\n",
    sum((is.na(res$L1) & res$nep1 & !is.na(res$L2)) |
        (is.na(res$L2) & res$nep2 & !is.na(res$L1)))))
cat(sprintf("  braids now dropped that were previously single-end labelled: %d\n",
    sum((is.na(res$L1) & !res$nep1 & !is.na(res$L2)) |
        (is.na(res$L2) & !res$nep2 & !is.na(res$L1)))))
cat("\n  regions parsed from blkID:", paste(head(unique(res$region),4), collapse=", "), "\n")
cat("  braids by class:\n"); print(table(res$cls))
BR$cls <- res$cls[match(BR$blkID, res$blkID)]
BR$cls[is.na(BR$cls)] <- "no call"

hr("4. DRAW")
lum <- function(h_,s_,v_){h <- rgb2hsv(col2rgb(h_)); hsv(h["h",],pmin(1,h["s",]*s_),pmin(1,h["v",]*v_))}
cols <- unique(BR$color)
Ac <- setNames(lum(cols,0.55,1.10), cols); Bc <- setNames(lum(cols,1.00,0.88), cols)
BR$fillcol <- ifelse(BR$cls=="A", Ac[BR$color], Bc[BR$color])
IV$barcol <- ifelse(IV$label=="A", "#1D9E75", "#D85A30")
PAL <- c(setNames(unique(BR$fillcol),unique(BR$fillcol)),
         setNames(unique(IV$barcol), unique(IV$barcol)))
p <- P
mkl <- function(d,al) layer(geom="polygon", stat="identity", position="identity",
  data=d, mapping=aes(x=x,y=y,group=blkID,fill=fillcol),
  params=list(alpha=al,colour=NA,linewidth=0), show.legend=FALSE)
p$layers <- c(list(mkl(BR[BR$cls=="A",],0.55), mkl(BR[BR$cls=="B",],0.88)),
              P$layers[-1],
              list(layer(geom="rect", stat="identity", position="identity",
                data=IV, mapping=aes(xmin=xa,xmax=xb,ymin=ybot,ymax=ytop,
                                     fill=barcol),
                params=list(colour=NA), show.legend=FALSE)))
p$scales$scales <- p$scales$scales[
  !vapply(p$scales$scales, function(sc) "fill" %in% sc$aesthetics, logical(1))]
p <- p + scale_fill_manual(values=PAL, guide="none") + theme_minimal(11) +
  theme(panel.background=element_rect(fill="grey96",colour=NA),
        plot.background=element_rect(fill="grey96",colour=NA),
        panel.grid=element_blank(), axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_text(face="italic",size=11),
        plot.title=element_text(face="bold",size=14),
        plot.subtitle=element_text(size=9,colour="grey30")) +
  labs(title="Droseraceae riparian, phased by subgenome",
       subtitle=paste0("Bars beneath each chromosome are the phased tracts: ",
         "green = A, orange = B. Braids are filtered against the same\n",
         "tracts -- drawn only where both ends fall on one subgenome. ",
         "MUTED = A, VIVID = B. Hue = ancestral region."),
       x="Chromosomes scaled by gene rank order", y=NULL)

hr("CHIMERIC BRAID ENDS — blocks that span an A/B transition")
cat("  A syntenic block can cross a subgenome junction: GENESPACE does not\n")
cat("  know about A/B, and homeologs are collinear. Taking a majority label\n")
cat("  hides that. This counts how often a braid end is genuinely mixed.\n\n")
comp <- function(x0, x1, y, reg) {
  k <- which(abs(ALLG$yrow - y) < 0.25 & ALLG$x >= x0 & ALLG$x <= x1 &
             (is.na(ALLG$reg) | ALLG$reg == reg))
  if (!length(k)) return(c(NA, NA, NA))
  t <- table(ALLG$lab[k])
  c(sum(t), max(t)/sum(t), length(t))
}
C1 <- t(mapply(comp, res$x1a, res$x1b, res$y1, res$region))
C2 <- t(mapply(comp, res$x2a, res$x2b, res$y2, res$region))
res$n1 <- C1[,1]; res$pur1 <- C1[,2]
res$n2 <- C2[,1]; res$pur2 <- C2[,2]
ok <- !is.na(res$pur1) & !is.na(res$pur2)
cat(sprintf("  braid ends assessable: %d of %d\n", sum(ok), nrow(res)))
for (thr in c(1.0, 0.99, 0.95, 0.90, 0.80)) {
  pure <- ok & res$pur1 >= thr & res$pur2 >= thr
  cat(sprintf("    both ends >= %.0f%% one label: %4d (%.1f%%)\n",
              100*thr, sum(pure), 100*mean(pure[ok])))
}
mixed <- ok & (res$pur1 < 0.95 | res$pur2 < 0.95)
cat(sprintf("\n  CHIMERIC (an end below 95%% pure): %d (%.1f%% of assessable)\n",
            sum(mixed), 100*mean(mixed[ok])))
if (sum(mixed)) {
  cat("  by region:\n")
  print(head(as.data.frame(sort(table(res$region[mixed]), decreasing=TRUE)), 8))
  cat("\n  A block spanning an A/B junction is REAL -- it means collinearity\n")
  cat("  continues across the transition. Drawing it as one colour is wrong.\n")
}
cat("\n  Current behaviour: majority label at 0.6 threshold in endlab().\n")
cat("  Braids below the purity cut are now DROPPED rather than coloured.\n")

hr("AUDIT — does every drawn braid match the bar beneath it?")
barlab_at <- function(x, y) {
  k <- which(abs(IV$y - y) < 0.25 & IV$xa <= x & IV$xb >= x)
  if (!length(k)) return(NA_character_)
  IV$label[k[1]]
}
D <- res[res$cls %in% c("A","B"), ]
D$bar1 <- mapply(barlab_at, (D$x1a+D$x1b)/2, D$y1)
D$bar2 <- mapply(barlab_at, (D$x2a+D$x2b)/2, D$y2)
D$mm1 <- !is.na(D$bar1) & D$bar1 != D$cls
D$mm2 <- !is.na(D$bar2) & D$bar2 != D$cls
cat(sprintf("  drawn braids: %d\n", nrow(D)))
cat(sprintf("  end 1 disagrees with its bar: %d (%.1f%%)\n",
            sum(D$mm1), 100*mean(D$mm1)))
cat(sprintf("  end 2 disagrees with its bar: %d (%.1f%%)\n",
            sum(D$mm2), 100*mean(D$mm2)))
MM <- D[D$mm1 | D$mm2, ]
if (nrow(MM)) {
  cat("\n  A mismatch is EXPECTED where the braid end spans two regions: the\n")
  cat("  bar shows the region at the midpoint, the braid reads its own region.\n")
  cat("  It is a BUG only if the braid's own region has no genes there.\n\n")
  cat("  worst offenders by region:\n")
  print(head(as.data.frame(sort(table(MM$region), decreasing=TRUE)), 8))
  write.csv(MM[, c("blkID","region","cls","bar1","bar2","y1","y2")],
            "DR/out/DR17_bar_braid_mismatch.csv", row.names=FALSE)
  cat("\n  wrote DR/out/DR17_bar_braid_mismatch.csv\n")
} else cat("\n  no mismatches -- bars and braids agree everywhere\n")

ggsave("DR/fig/DR17_riparian_AB.pdf", p, width=16, height=9)
ggsave("DR/fig/DR17_riparian_AB.png", p, width=16, height=9, dpi=200)
cat("  wrote DR/fig/DR17_riparian_AB.{pdf,png}\n")
cat("\n  Read section 3 first: a large 'no call' count means the figure shows\n")
cat("  missing data as much as biology.\n")
