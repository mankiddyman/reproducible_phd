#!/usr/bin/env Rscript
# ============================================================================
# DR22 — the riparian in PHYSICAL coordinates, + the occlusion test
#
# Uses Nepenthes_gracilis_bp_rSourceData.rda (same GENESPACE run, Aug 11).
# Section 1 gates: bp -> x must be affine per chromosome, or we stop.
# Section 5 tests whether the visual defect is BRAID CROSSING, not placement:
#   bars occupy y-0.05..y-0.11, exactly where a descending ribbon passes.
# Two bar styles drawn so the cause is visible rather than argued.
#
# OUT DR/fig/DR22_bp_band.pdf, DR/fig/DR22_bp_boxfill.pdf,
#     DR/out/DR22_occlusion.csv
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({ library(GENESPACE); library(data.table); library(ggplot2) })
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=200)
GS <- file.path(dirname(getwd()), "genespace")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }

hr("1. LOAD bp OBJECT AND FIT bp -> PLOT X")
e <- new.env(); load(file.path(GS,"riparian","Nepenthes_gracilis_bp_rSourceData.rda"), envir=e)
P <- e$srcd$ggplotObj
BR    <- as.data.frame(P$layers[[1]]$data)
CHROM <- as.data.frame(e$srcd$sourceData$chromosomes); CHROM$key <- paste(CHROM$genome, CHROM$chr)
BLK   <- as.data.frame(e$srcd$sourceData$blocks)
cat(sprintf("  braid rows %d | blocks %d | chromosomes %d\n", nrow(BR), nrow(BLK), nrow(CHROM)))
cat("  chromosomes cols:", paste(names(CHROM), collapse=", "), "\n")
cat(sprintf("  box width vs length: cor %.6f | ratio range %.6g .. %.6g\n",
    cor(CHROM$x2-CHROM$x1, CHROM$length),
    min((CHROM$x2-CHROM$x1)/CHROM$length), max((CHROM$x2-CHROM$x1)/CHROM$length)))

EP <- as.data.frame(as.data.table(BR)[, .(xlo=min(x), xhi=max(x)), by=.(blkID, yrow=round(y,6))])
SD <- rbind(
  data.frame(blkID=BLK$blkID, yrow=round(as.numeric(BLK$y1),6), key=paste(BLK$genome1,BLK$chr1),
             bpa=BLK$startBp1, bpb=BLK$endBp1, stringsAsFactors=FALSE),
  data.frame(blkID=BLK$blkID, yrow=round(as.numeric(BLK$y2),6), key=paste(BLK$genome2,BLK$chr2),
             bpa=BLK$startBp2, bpb=BLK$endBp2, stringsAsFactors=FALSE))
SD <- merge(SD, EP, by=c("blkID","yrow"))
cat(sprintf("  block sides with a polygon attachment: %d\n", nrow(SD)))

FIT <- do.call(rbind, lapply(split(SD, SD$key), function(d){
  b <- c(pmin(d$bpa,d$bpb), pmax(d$bpa,d$bpb))
  if (length(unique(b)) < 2) return(data.frame(key=d$key[1], n=nrow(d), r2=NA, a=NA, s=NA, mx=NA))
  f <- function(x){ m <- lm(x ~ b)
    c(summary(m)$r.squared, unname(coef(m)[1]), unname(coef(m)[2]), max(abs(residuals(m)))) }
  p <- f(c(d$xlo,d$xhi)); q <- f(c(d$xhi,d$xlo)); w <- if (p[1] >= q[1]) p else q
  data.frame(key=d$key[1], n=nrow(d), r2=w[1], a=w[2], s=w[3], mx=w[4]) }))
FIT$mb_err <- FIT$mx / abs(FIT$s) / 1e6
cat(sprintf("\n  worst R^2 %.12g | chromosomes R^2<0.9999: %d | slopes negative: %d\n",
    min(FIT$r2, na.rm=TRUE), sum(FIT$r2 < 0.9999, na.rm=TRUE), sum(FIT$s < 0, na.rm=TRUE)))
cat("  residual expressed in Mb:\n"); print(summary(FIT$mb_err))
print(head(FIT[order(-FIT$mb_err), c("key","n","r2","s","mx","mb_err")], 8), row.names=FALSE)
if (any(is.na(FIT$r2))) cat(sprintf("  chromosomes with too few blocks to fit: %d\n", sum(is.na(FIT$r2))))
if (max(FIT$mb_err, na.rm=TRUE) > 0.05) die("bp -> x not affine; worst residual %.3f Mb", max(FIT$mb_err, na.rm=TRUE))
cat("  bp -> x is affine. Proceeding.\n")

hr("2. GENES -> PLOT X, LABELS")
BED <- as.data.frame(fread(file.path(GS,"results","combBed.txt")))
BED$key <- paste(BED$genome, BED$chr); BED <- BED[BED$key %in% CHROM$key, ]
BED$bp  <- (BED$start + BED$end)/2
BED$px  <- FIT$a[match(BED$key, FIT$key)] + FIT$s[match(BED$key, FIT$key)] * BED$bp
BED$cx1 <- pmin(CHROM$x1, CHROM$x2)[match(BED$key, CHROM$key)]
BED$cx2 <- pmax(CHROM$x1, CHROM$x2)[match(BED$key, CHROM$key)]
BED <- BED[!is.na(BED$px), ]
BED$px <- pmin(pmax(BED$px, BED$cx1), BED$cx2)
G <- read.csv("DR/out/DR02_propagated_blocks.csv", stringsAsFactors=FALSE)
G$key <- paste(G$genome, G$chr); G <- G[G$label %in% c("A","B"), ]
m <- match(paste(BED$key,BED$id), paste(G$key,G$gene))
BED$lab <- G$label[m]; BED$reg <- G$region[m]
PAIRS <- read.csv("fractionation_by_chrpair.csv", stringsAsFactors=FALSE)
s1 <- PAIRS$retained_more == PAIRS$chrA
DIOSIDE <- c(setNames(ifelse(s1,"A","B"), PAIRS$chrA), setNames(ifelse(s1,"B","A"), PAIRS$chrB))
di <- BED$genome == "Dionaea_muscipula"
BED$lab[di] <- unname(DIOSIDE[BED$chr[di]]); BED$reg[di] <- "Dionaea"
cat(sprintf("  genes placed %d | labelled %d\n", nrow(BED), sum(!is.na(BED$lab))))

hr("3. TRACTS")
CHR <- data.frame(genome=CHROM$genome, chr=CHROM$chr, key=CHROM$key,
                  x0=pmin(CHROM$x1,CHROM$x2), x1=pmax(CHROM$x1,CHROM$x2),
                  cy1=pmin(CHROM$y1,CHROM$y2), cy2=pmax(CHROM$y1,CHROM$y2),
                  y=(CHROM$y1+CHROM$y2)/2, stringsAsFactors=FALSE)
LG <- BED[!is.na(BED$lab) & !di, ]
TR <- do.call(rbind, lapply(split(LG, LG$key), function(z){
  z <- z[order(z$px), ]; r <- rle(paste(z$reg, z$lab))
  en <- cumsum(r$lengths); st <- en - r$lengths + 1
  do.call(rbind, lapply(seq_along(r$lengths), function(j){ w <- z[st[j]:en[j], ]
    data.frame(genome=w$genome[1], chr=w$chr[1], region=w$reg[1], label=w$lab[1],
               xa=min(w$px), xb=max(w$px), n=nrow(w), stringsAsFactors=FALSE) })) }))
TR <- TR[TR$n >= 10, ]
.r5 <- TR[TR$genome=="Drosera_regia" & TR$chr=="chr5_collapsed", ]
if (length(unique(.r5$label)) != 2) { print(.r5); die("regia chr5_collapsed lost its A/B split") }
.both <- sum(vapply(split(TR, paste(TR$genome,TR$chr)), function(z) length(unique(z$label))>1, logical(1)))
if (.both < 25) die("only %d chromosomes carry both labels", .both)
TR <- merge(TR, CHR[, c("genome","chr","y","cy1","cy2")], by=c("genome","chr"))
DC <- CHR[CHR$genome=="Dionaea_muscipula", ]; DC$label <- unname(DIOSIDE[DC$chr])
IV <- rbind(TR[, c("genome","chr","region","label","xa","xb","y","cy1","cy2")],
            data.frame(genome=DC$genome, chr=DC$chr, region="Dionaea", label=DC$label,
                       xa=DC$x0, xb=DC$x1, y=DC$y, cy1=DC$cy1, cy2=DC$cy2, stringsAsFactors=FALSE))
cat(sprintf("  tracts %d (%d A, %d B) | both-label chromosomes %d | bars %d\n",
    nrow(TR), sum(TR$label=="A"), sum(TR$label=="B"), .both, nrow(IV)))

hr("4. BRAID CLASS + ATTACHMENT AUDIT (span overlap, region-agnostic)")
SD$L <- NA_character_
for (i in seq_len(nrow(SD))) {
  z <- BED$lab[BED$key==SD$key[i] & BED$bp >= min(SD$bpa[i],SD$bpb[i]) & BED$bp <= max(SD$bpa[i],SD$bpb[i])]
  z <- z[!is.na(z)]; if (!length(z)) next
  pu <- max(sum(z=="A"), sum(z=="B"))/length(z)
  if (pu >= 0.95) SD$L[i] <- if (sum(z=="A") >= sum(z=="B")) "A" else "B" }
W <- merge(SD[!duplicated(paste(SD$blkID,SD$yrow)), ], SD, by="blkID", suffixes=c("1","2"))
W <- W[W$yrow1 < W$yrow2, ]; W <- W[!duplicated(W$blkID), ]
nep <- function(k) startsWith(k, "Nepenthes_gracilis ")
W$cls <- ifelse(!is.na(W$L1) & !is.na(W$L2), ifelse(W$L1==W$L2, W$L1, "discordant"),
          ifelse(is.na(W$L1) & nep(W$key1) & !is.na(W$L2), W$L2,
          ifelse(is.na(W$L2) & nep(W$key2) & !is.na(W$L1), W$L1, "no call")))
cat("  braids by class:\n"); print(table(W$cls))
ovl <- function(x0,x1,y,cls){ k <- which(abs(IV$y-y)<0.25 & IV$xb>x0 & IV$xa<x1)
  if(!length(k)) return(0); w <- pmin(IV$xb[k],x1)-pmax(IV$xa[k],x0)
  sum(w[IV$label[k] != cls])/max(x1-x0,1e-9) }
D <- W[W$cls %in% c("A","B"), ]
D$f1 <- mapply(ovl, D$xlo1, D$xhi1, D$yrow1, D$cls)
D$f2 <- mapply(ovl, D$xlo2, D$xhi2, D$yrow2, D$cls)
D$fmax <- pmax(D$f1, D$f2)
cat(sprintf("  candidates %d | ends touching a wrong-coloured bar: %d (%.1f%%)\n",
    nrow(D), sum(D$fmax > 0), 100*mean(D$fmax > 0)))
KEEP <- D$blkID[D$fmax <= 0]
cat(sprintf("  kept %d | dropped %d\n", length(KEEP), nrow(D)-length(KEEP)))

hr("5. OCCLUSION TEST — where does the ribbon sit AT THE BAR BAND?")
cat("  Bars occupy y-0.05..y-0.11. A ribbon descending from that row passes\n")
cat("  straight through. If braids cross wrong-coloured bars there, the defect\n")
cat("  is CROSSING, not placement -- and no coordinate system fixes it.\n\n")
BRdt <- as.data.table(BR)
occ <- rbindlist(lapply(KEEP, function(id){
  v <- BRdt[blkID==id]; ytop <- max(v$y)
  s <- v[y <= ytop-0.045 & y >= ytop-0.125]
  if (!nrow(s)) return(NULL)
  cl <- D$cls[match(id, D$blkID)]
  data.table(blkID=id, cls=cl, yrow=ytop, bx0=min(s$x), bx1=max(s$x),
             att0=min(v$x[abs(v$y-ytop)<1e-9]), att1=max(v$x[abs(v$y-ytop)<1e-9])) }))
occ[, cross := mapply(ovl, bx0, bx1, yrow, cls)]
occ[, att   := mapply(ovl, att0, att1, yrow, cls)]
occ[, drift := pmin(abs(bx0-att0), abs(bx1-att1))]
cat(sprintf("  upper attachments tested: %d\n", nrow(occ)))
cat("  fraction of the ribbon's width over a WRONG-coloured bar, at bar depth:\n")
print(summary(occ$cross))
for (th in c(0, 0.05, 0.25, 0.5)) cat(sprintf("    exceeding %3.0f%%: %4d (%.1f%%)\n",
  100*th, sum(occ$cross > th), 100*mean(occ$cross > th)))
cat(sprintf("\n  by comparison, AT the attachment line itself: %d of %d (%.1f%%) wrong\n",
    sum(occ$att > 0), nrow(occ), 100*mean(occ$att > 0)))
cat("  horizontal drift between attachment and bar depth (plot units):\n"); print(summary(occ$drift))
fwrite(occ, "DR/out/DR22_occlusion.csv")
cat("\n  VERDICT: ")
if (mean(occ$cross > 0) > 5*max(mean(occ$att > 0), 0.001))
  cat("CROSSING CONFIRMED -- braids are correct where they land and wrong\n           where they cross. Use the box-fill figure.\n") else
  cat("crossing is NOT the dominant effect; look elsewhere.\n")

hr("6. DRAW — two bar styles")
lum <- function(h_,s_,v_){h <- rgb2hsv(col2rgb(h_)); hsv(h["h",],pmin(1,h["s",]*s_),pmin(1,h["v",]*v_))}
BR$cls <- W$cls[match(BR$blkID, W$blkID)]
BR$cls[is.na(BR$cls) | !(BR$blkID %in% KEEP)] <- "no call"
cols <- unique(BR$color)
Ac <- setNames(lum(cols,0.55,1.10), cols); Bc <- setNames(lum(cols,1.00,0.88), cols)
BR$fillcol <- ifelse(BR$cls=="A", Ac[BR$color], Bc[BR$color])
IV$barcol <- ifelse(IV$label=="A", "#1D9E75", "#D85A30")
PAL <- c(setNames(unique(BR$fillcol),unique(BR$fillcol)), setNames(unique(IV$barcol),unique(IV$barcol)))
mkl <- function(d,al) layer(geom="polygon", stat="identity", position="identity", data=d,
  mapping=aes(x=x,y=y,group=blkID,fill=fillcol), params=list(alpha=al,colour=NA,linewidth=0), show.legend=FALSE)
base <- function(sub) theme_minimal(11) + theme(
  panel.background=element_rect(fill="grey96",colour=NA), plot.background=element_rect(fill="grey96",colour=NA),
  panel.grid=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(),
  axis.text.y=element_text(face="italic",size=11), plot.title=element_text(face="bold",size=14),
  plot.subtitle=element_text(size=9,colour="grey30"))
lab <- function(sub) labs(title="Droseraceae riparian, phased by subgenome",
  subtitle=sub, x="Chromosomes scaled by PHYSICAL LENGTH (bp)", y=NULL)

IVb <- IV; IVb$ytop <- IVb$y-0.05; IVb$ybot <- IVb$y-0.11
p1 <- P; p1$layers <- c(list(mkl(BR[BR$cls=="A",],0.55), mkl(BR[BR$cls=="B",],0.88)), P$layers[-1],
  list(layer(geom="rect", stat="identity", position="identity", data=IVb,
    mapping=aes(xmin=xa,xmax=xb,ymin=ybot,ymax=ytop,fill=barcol), params=list(colour=NA), show.legend=FALSE)))
p1$scales$scales <- p1$scales$scales[!vapply(p1$scales$scales, function(s) "fill" %in% s$aesthetics, logical(1))]
p1 <- p1 + scale_fill_manual(values=PAL, guide="none") + base() +
  lab("Bars BENEATH each chromosome (green = A, orange = B). Physical scale.")
ggsave("DR/fig/DR22_bp_band.pdf", p1, width=16, height=9)

h <- IV$cy2 - IV$cy1
IVf <- IV; IVf$ybot <- IV$cy1 + 0.18*h; IVf$ytop <- IV$cy2 - 0.18*h
p2 <- P; p2$layers <- c(list(mkl(BR[BR$cls=="A",],0.55), mkl(BR[BR$cls=="B",],0.88)), P$layers[-1],
  list(layer(geom="rect", stat="identity", position="identity", data=IVf,
    mapping=aes(xmin=xa,xmax=xb,ymin=ybot,ymax=ytop,fill=barcol), params=list(colour=NA), show.legend=FALSE)))
p2$scales$scales <- p2$scales$scales[!vapply(p2$scales$scales, function(s) "fill" %in% s$aesthetics, logical(1))]
p2 <- p2 + scale_fill_manual(values=PAL, guide="none") + base() +
  lab("Tracts filled INTO each chromosome (green = A, orange = B), so no braid can cross them. Physical scale.")
ggsave("DR/fig/DR22_bp_boxfill.pdf", p2, width=16, height=9)
cat(sprintf("  wrote DR/fig/DR22_bp_band.pdf and DR/fig/DR22_bp_boxfill.pdf | %d braids drawn\n", length(KEEP)))
