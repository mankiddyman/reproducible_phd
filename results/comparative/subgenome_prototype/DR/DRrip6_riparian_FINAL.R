#!/usr/bin/env Rscript
# ============================================================================
# DR25 — PRODUCTION riparian: gene rank, edge strips, VERTICAL ATTACHMENT STEMS
#
# Coordinate map (DR18/19/20, exact to 0.000 plot units over 2556 block sides):
#     plot_ord = rank(ord) within chromosome - offset(chromosome)
#     x        = chromosome$x1 + plot_ord          [no chromosome is flipped]
#
# DR24 fixed WHERE braids attach. This fixes what the attachment LOOKS like.
# A ribbon left the box at a shallow angle, so a few pixels below the edge it
# was already under a distant part of the chromosome -- e.g. capensis chr15:
# blue chr6_dom/B braids attach at 64.6-84.8% (orange) but immediately pass
# beneath 37.2-56.3% (green chr7_dom/A). The eye reads the crossing as the
# landing.
#
# Each polygon's y-range is compressed by STEM at both ends (monotone remap --
# ring order and curve shape preserved), and a rectangle of the same fill and
# alpha spans the attachment's exact x-range over that STEM. The ribbon now
# leaves vertically, so pixels beneath a strip belong only to braids that
# attach there.
#
# OUT DR/fig/DR25_riparian_AB.{pdf,png}
#     DR/fig/DR25_zoom_{capensis_chr15,regia_chr5,scorpioides_chr1}.png
#     DR/out/DR25_dropped.csv, DR/out/DR25_tracts.csv
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({ library(GENESPACE); library(data.table); library(ggplot2) })
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
GS <- file.path(dirname(getwd()), "genespace")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }

BARA <- "#1D9E75"    # A subgenome strip
BARB <- "#D85A30"    # B subgenome strip
FR   <- 0.30         # strip height, fraction of chromosome box
HAIR <- 0.012        # white hairline between strip and braid
STEM <- 0.12         # vertical riser at each attachment  <-- the knob

hr("1. COORDINATE MAP")
e <- new.env(); load(file.path(GS,"riparian","Nepenthes_gracilis_geneOrder_rSourceData.rda"), envir=e)
P <- e$srcd$ggplotObj
BR    <- as.data.frame(P$layers[[1]]$data)
CHROM <- as.data.frame(e$srcd$sourceData$chromosomes); CHROM$key <- paste(CHROM$genome, CHROM$chr)
BLK   <- as.data.frame(e$srcd$sourceData$blocks)
cat(sprintf("  saved plot layers %d | braid rows %d | blocks %d\n", length(P$layers), nrow(BR), nrow(BLK)))
BED <- as.data.frame(fread(file.path(GS,"results","combBed.txt")))
BED$key <- paste(BED$genome, BED$chr); BED <- BED[BED$key %in% CHROM$key, ]
BED <- BED[order(BED$key, BED$ord), ]; BED$r <- ave(BED$ord, BED$key, FUN=function(z) z-min(z)+1L)
AN <- unique(rbind(
  data.frame(key=paste(BLK$genome1,BLK$chr1), gene=BLK$firstGene1, pord=BLK$startOrd1),
  data.frame(key=paste(BLK$genome1,BLK$chr1), gene=BLK$lastGene1,  pord=BLK$endOrd1),
  data.frame(key=paste(BLK$genome2,BLK$chr2), gene=BLK$firstGene2, pord=BLK$startOrd2),
  data.frame(key=paste(BLK$genome2,BLK$chr2), gene=BLK$lastGene2,  pord=BLK$endOrd2), stringsAsFactors=FALSE))
AN$r <- BED$r[match(paste(AN$key,AN$gene), paste(BED$key,BED$ofID))]; AN <- AN[!is.na(AN$r), ]
if (max(tapply(AN$r-AN$pord, AN$key, function(z) diff(range(z)))) != 0) die("offset not constant")
OFF <- tapply(AN$r-AN$pord, AN$key, function(z) z[1])
BED$cx1 <- CHROM$x1[match(BED$key,CHROM$key)]; BED$clen <- CHROM$length[match(BED$key,CHROM$key)]
BED$px  <- BED$cx1 + pmin(pmax(BED$r - OFF[BED$key], 1), BED$clen)
EP <- as.data.frame(as.data.table(BR)[, .(xlo=min(x), xhi=max(x)), by=.(blkID, yrow=round(y,6))])
SD <- rbind(
  data.frame(blkID=BLK$blkID, yrow=round(as.numeric(BLK$y1),6), key=paste(BLK$genome1,BLK$chr1),
             g1=BLK$firstGene1, g2=BLK$lastGene1, stringsAsFactors=FALSE),
  data.frame(blkID=BLK$blkID, yrow=round(as.numeric(BLK$y2),6), key=paste(BLK$genome2,BLK$chr2),
             g1=BLK$firstGene2, g2=BLK$lastGene2, stringsAsFactors=FALSE))
SD <- merge(SD, EP, by=c("blkID","yrow"))
kk <- paste(BED$key, BED$ofID)
SD$xa <- BED$px[match(paste(SD$key,SD$g1),kk)]; SD$xb <- BED$px[match(paste(SD$key,SD$g2),kk)]
err <- max(pmax(abs(pmin(SD$xa,SD$xb)-SD$xlo), abs(pmax(SD$xa,SD$xb)-SD$xhi)), na.rm=TRUE)
cat(sprintf("  block sides %d | max coordinate error %.4g plot units\n", nrow(SD), err))
if (err > 1e-6) die("coordinate map regressed")

hr("2. GENE LABELS")
G <- read.csv("DR/out/DR02_propagated_blocks.csv", stringsAsFactors=FALSE)
G$key <- paste(G$genome, G$chr); G <- G[G$label %in% c("A","B"), ]
m <- match(paste(BED$key,BED$id), paste(G$key,G$gene)); BED$lab <- G$label[m]; BED$reg <- G$region[m]
PAIRS <- read.csv("fractionation_by_chrpair.csv", stringsAsFactors=FALSE); s1 <- PAIRS$retained_more==PAIRS$chrA
DIOSIDE <- c(setNames(ifelse(s1,"A","B"),PAIRS$chrA), setNames(ifelse(s1,"B","A"),PAIRS$chrB))
di <- BED$genome=="Dionaea_muscipula"
BED$lab[di] <- unname(DIOSIDE[BED$chr[di]]); BED$reg[di] <- "Dionaea"
cat(sprintf("  labelled genes %d of %d\n", sum(!is.na(BED$lab)), nrow(BED)))

hr("3. TRACTS")
CHR <- data.frame(key=CHROM$key, genome=CHROM$genome, chr=CHROM$chr,
                  x0=CHROM$x1, x1=CHROM$x2, by1=pmin(CHROM$y1,CHROM$y2),
                  by2=pmax(CHROM$y1,CHROM$y2), y=(CHROM$y1+CHROM$y2)/2, stringsAsFactors=FALSE)
LG <- BED[!is.na(BED$lab) & !di, ]
TR <- do.call(rbind, lapply(split(LG, LG$key), function(z){ z <- z[order(z$px), ]
  r <- rle(paste(z$reg,z$lab)); en <- cumsum(r$lengths); st <- en-r$lengths+1
  do.call(rbind, lapply(seq_along(r$lengths), function(j){ w <- z[st[j]:en[j], ]
    data.frame(key=w$key[1], region=w$reg[1], label=w$lab[1],
               xa=min(w$px), xb=max(w$px), n=nrow(w), stringsAsFactors=FALSE) })) }))
TR <- TR[TR$n >= 10, ]
TRk <- merge(TR, CHR, by="key")
.r5 <- TRk[TRk$key=="Drosera_regia chr5_collapsed", ]
if (length(unique(.r5$label)) != 2) { print(.r5); die("regia chr5_collapsed lost its A/B split") }
.both <- sum(vapply(split(TRk, TRk$key), function(z) length(unique(z$label))>1, logical(1)))
if (.both < 25) die("only %d chromosomes carry both labels", .both)
DC <- CHR[CHR$genome=="Dionaea_muscipula", ]
IV <- rbind(TRk[, c("key","genome","chr","region","label","xa","xb","n","x0","x1","by1","by2","y")],
            data.frame(key=DC$key, genome=DC$genome, chr=DC$chr, region="Dionaea",
                       label=unname(DIOSIDE[DC$chr]), xa=DC$x0, xb=DC$x1, n=NA,
                       x0=DC$x0, x1=DC$x1, by1=DC$by1, by2=DC$by2, y=DC$y, stringsAsFactors=FALSE))
cat(sprintf("  tracts %d (%d A, %d B) | both-label chromosomes %d | strips %d\n",
            nrow(TR), sum(TR$label=="A"), sum(TR$label=="B"), .both, nrow(IV)))
write.csv(IV, "DR/out/DR25_tracts.csv", row.names=FALSE)

hr("4. BRAID CLASSES FROM EACH BLOCK'S OWN GENES")
SD$ra <- BED$r[match(paste(SD$key,SD$g1),kk)]; SD$rb <- BED$r[match(paste(SD$key,SD$g2),kk)]
DT <- as.data.table(BED[,c("key","r","lab")]); setkey(DT,key,r)
SD$L <- NA_character_
for (i in seq_len(nrow(SD))) { if (is.na(SD$ra[i])||is.na(SD$rb[i])) next
  z <- DT[.(SD$key[i])][r>=min(SD$ra[i],SD$rb[i]) & r<=max(SD$ra[i],SD$rb[i]), lab]; z <- z[!is.na(z)]
  if (!length(z)) next
  if (max(sum(z=="A"),sum(z=="B"))/length(z) >= 0.95) SD$L[i] <- if (sum(z=="A")>=sum(z=="B")) "A" else "B" }
W <- merge(SD, SD, by="blkID", suffixes=c("1","2")); W <- W[W$yrow1 < W$yrow2, ]; W <- W[!duplicated(W$blkID),]
nep <- function(k) startsWith(k,"Nepenthes_gracilis ")
W$cls <- ifelse(!is.na(W$L1)&!is.na(W$L2), ifelse(W$L1==W$L2,W$L1,"discordant"),
          ifelse(is.na(W$L1)&nep(W$key1)&!is.na(W$L2), W$L2,
          ifelse(is.na(W$L2)&nep(W$key2)&!is.na(W$L1), W$L1, "no call")))
cat("  braids by class:\n"); print(table(W$cls))

hr("5. ENFORCE — attachment must sit wholly on its own colour")
ovl <- function(key,x0,x1,cls){ k <- which(IV$key==key & IV$xb>x0 & IV$xa<x1)
  if(!length(k)) return(0); w <- pmin(IV$xb[k],x1)-pmax(IV$xa[k],x0); sum(w[IV$label[k]!=cls])/max(x1-x0,1e-9) }
D <- W[W$cls %in% c("A","B"), ]
D$f1 <- mapply(ovl, D$key1, D$xlo1, D$xhi1, D$cls)
D$f2 <- mapply(ovl, D$key2, D$xlo2, D$xhi2, D$cls)
D$fmax <- pmax(D$f1,D$f2)
KEEP <- D$blkID[D$fmax <= 0]; DROP <- D[D$fmax > 0, ]
cat(sprintf("  candidates %d | kept %d | dropped %d (%.1f%%)\n",
            nrow(D), length(KEEP), nrow(DROP), 100*mean(D$fmax>0)))
if (nrow(DROP)) write.csv(DROP[,c("blkID","key1","key2","cls","f1","f2","fmax")],
                          "DR/out/DR25_dropped.csv", row.names=FALSE)
BR$cls <- W$cls[match(BR$blkID,W$blkID)]; BR$cls[is.na(BR$cls) | !(BR$blkID %in% KEEP)] <- "no call"
if (nrow(D[D$blkID %in% KEEP, ]) && max(D$fmax[D$blkID %in% KEEP]) > 0) die("enforcement failed")
cat("  GUARANTEE: every drawn attachment lies wholly on strips of its own subgenome.\n")

hr("6. STEMS — compress each ribbon, abut a vertical riser at each end")
lum <- function(h_,s_,v_){h <- rgb2hsv(col2rgb(h_)); hsv(h["h",],pmin(1,h["s",]*s_),pmin(1,h["v",]*v_))}
cols <- unique(BR$color); Ac <- setNames(lum(cols,0.55,1.10),cols); Bc <- setNames(lum(cols,1.00,0.88),cols)
BR$fillcol <- ifelse(BR$cls=="A", Ac[BR$color], Bc[BR$color])
Bd <- as.data.table(BR[BR$cls %in% c("A","B"), ])
Bd[, `:=`(ylo = min(y), yhi = max(y)), by = blkID]
Bd[, span := yhi - ylo]
cat("  ribbon vertical spans:\n"); print(summary(unique(Bd[, .(blkID, span)])$span))
if (min(Bd$span) <= 2*STEM) die("STEM %.3f too large for the shortest ribbon (%.3f)", STEM, min(Bd$span))
STM <- Bd[, .(xlo_b = min(x[abs(y-ylo) < 1e-9]), xhi_b = max(x[abs(y-ylo) < 1e-9]),
              xlo_t = min(x[abs(y-yhi) < 1e-9]), xhi_t = max(x[abs(y-yhi) < 1e-9]),
              ylo = ylo[1], yhi = yhi[1]), by = .(blkID, cls, fillcol)]
Bd[, y := ylo + STEM + (y - ylo) * (span - 2*STEM)/span]        # monotone remap
BRp <- as.data.frame(Bd)
S1 <- data.frame(blkID=STM$blkID, cls=STM$cls, fillcol=STM$fillcol,
                 xa=STM$xlo_b, xb=STM$xhi_b, ymin=STM$ylo, ymax=STM$ylo+STEM, stringsAsFactors=FALSE)
S2 <- data.frame(blkID=STM$blkID, cls=STM$cls, fillcol=STM$fillcol,
                 xa=STM$xlo_t, xb=STM$xhi_t, ymin=STM$yhi-STEM, ymax=STM$yhi, stringsAsFactors=FALSE)
SS <- rbind(S1, S2)
SS <- SS[SS$xb > SS$xa, ]
cat(sprintf("  stems built: %d (2 per drawn braid) | stem height %.3f | widths: ", nrow(SS), STEM))
cat(sprintf("median %.0f, max %.0f plot units\n", median(SS$xb-SS$xa), max(SS$xb-SS$xa)))
# every stem carries the key of the chromosome it stands on: lower stem (S1)
# is the ribbon's ylo end, upper stem (S2) the yhi end. W$yrow1 < W$yrow2, so
# key1 is the lower row and key2 the upper.
SS$key <- c(W$key1[match(S1$blkID, W$blkID)], W$key2[match(S2$blkID, W$blkID)])
SS <- SS[!is.na(SS$key), ]
SS$bad <- mapply(ovl, SS$key, SS$xa, SS$xb, SS$cls)
cat(sprintf("  stem footprints overlapping a wrong-coloured strip: %d of %d\n",
            sum(SS$bad > 0), nrow(SS)))
if (any(SS$bad > 0)) {
  b <- SS[SS$bad > 0, ]; b$width <- b$xb - b$xa
  cat("  offenders (stem rectangle wider than the polygon's own footprint):\n")
  print(head(b[order(-b$bad), c("key","cls","xa","xb","width","bad")], 12), row.names=FALSE)
  # clamp each stem to the widest wrong-free window inside its own footprint
  shrink <- function(key, x0, x1, cls) {
    k <- which(IV$key==key & IV$label==cls & IV$xb>x0 & IV$xa<x1)
    if (!length(k)) return(c(NA_real_, NA_real_))
    w <- pmin(IV$xb[k],x1) - pmax(IV$xa[k],x0); j <- k[which.max(w)]
    c(max(IV$xa[j], x0), min(IV$xb[j], x1)) }
  fix <- t(mapply(shrink, b$key, b$xa, b$xb, b$cls))
  SS$xa[SS$bad > 0] <- fix[,1]; SS$xb[SS$bad > 0] <- fix[,2]
  SS <- SS[!is.na(SS$xa) & SS$xb > SS$xa, ]
  SS$bad <- mapply(ovl, SS$key, SS$xa, SS$xb, SS$cls)
  cat(sprintf("  after clamping to the largest same-colour window: %d bad, %d stems remain\n",
              sum(SS$bad > 0), nrow(SS)))
}
if (any(SS$bad > 0)) die("%d stems still overlap a wrong-coloured strip", sum(SS$bad > 0))
cat("  GUARANTEE: every stem rectangle lies wholly on strips of its own subgenome.\n")

hr("7. DRAW")
h <- IV$by2 - IV$by1
ST <- rbind(transform(IV, ymin = by2 - FR*h, ymax = by2),
            transform(IV, ymin = by1,        ymax = by1 + FR*h))
ST$barcol <- ifelse(ST$label=="A", BARA, BARB)
SEP <- rbind(transform(IV, ymin=by2-HAIR, ymax=by2, xa=x0, xb=x1),
             transform(IV, ymin=by1, ymax=by1+HAIR, xa=x0, xb=x1))
SEP <- SEP[!duplicated(paste(SEP$key, SEP$ymin)), ]; SEP$barcol <- "#FFFFFF"
PAL <- c(setNames(unique(BRp$fillcol),unique(BRp$fillcol)),
         setNames(c(BARA,BARB,"#FFFFFF"), c(BARA,BARB,"#FFFFFF")))
poly <- function(d,al) layer(geom="polygon", stat="identity", position="identity", data=d,
  mapping=aes(x=x,y=y,group=blkID,fill=fillcol), params=list(alpha=al,colour=NA,linewidth=0), show.legend=FALSE)
srec <- function(d,al) layer(geom="rect", stat="identity", position="identity", data=d,
  mapping=aes(xmin=xa,xmax=xb,ymin=ymin,ymax=ymax,fill=fillcol), params=list(alpha=al,colour=NA), show.legend=FALSE)
brec <- function(d) layer(geom="rect", stat="identity", position="identity", data=d,
  mapping=aes(xmin=xa,xmax=xb,ymin=ymin,ymax=ymax,fill=barcol), params=list(colour=NA), show.legend=FALSE)
nL <- length(P$layers)
build <- function(bb, ss) { q <- P
  q$layers <- c(list(poly(bb[bb$cls=="A",],0.55), srec(ss[ss$cls=="A",],0.55),
                     poly(bb[bb$cls=="B",],0.88), srec(ss[ss$cls=="B",],0.88)),
                P$layers[2], list(brec(ST), brec(SEP)), if (nL>=3) P$layers[3:nL] else NULL)
  q$scales$scales <- q$scales$scales[!vapply(q$scales$scales,function(s)"fill"%in%s$aesthetics,logical(1))]
  q + scale_fill_manual(values=PAL, guide="none") }
sty <- function(q, ttl, sub) q + theme_minimal(11) + theme(
  panel.background=element_rect(fill="grey96",colour=NA), plot.background=element_rect(fill="grey96",colour=NA),
  panel.grid=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(),
  axis.text.y=element_text(face="italic",size=11), plot.title=element_text(face="bold",size=14),
  plot.subtitle=element_text(size=9,colour="grey30")) +
  labs(title=ttl, subtitle=sub, x="Chromosomes scaled by GENE RANK -- physical spacing is not preserved", y=NULL)
SUB <- paste0("Coloured strips along each chromosome edge are the phased subgenome tracts: green = A, orange = B. ",
  "Each ribbon leaves its\nchromosome vertically, so the colour it touches is the colour it attaches to. Drawn only ",
  "where both ends are >=95%\none subgenome among the block's own genes. MUTED = A, VIVID = B; hue = ancestral region, ",
  "independent of green/orange.")
p <- sty(build(BRp, SS), "Droseraceae riparian, phased by subgenome", SUB)
ggsave("DR/fig/DR25_riparian_AB.pdf", p, width=16, height=9)
ggsave("DR/fig/DR25_riparian_AB.png", p, width=16, height=9, dpi=200)
cat(sprintf("  wrote DR/fig/DR25_riparian_AB.{pdf,png} | %d braids drawn\n", length(KEEP)))

hr("8. ZOOM VERIFICATION")
zoom <- function(key, file){ i <- which(CHR$key==key); if(!length(i)) return(invisible())
  XL <- c(CHR$x0[i]-60, CHR$x1[i]+60); YL <- c(CHR$y[i]-0.55, CHR$y[i]+0.55)
  kb <- unique(BRp$blkID[BRp$x>=XL[1] & BRp$x<=XL[2] & abs(BRp$y-CHR$y[i])<0.55])
  q <- build(BRp[BRp$blkID %in% kb, ], SS[SS$blkID %in% kb, ])
  q <- sty(q, key, "zoom -- each ribbon leaves vertically; the colour it touches is where it attaches") +
       coord_cartesian(xlim=XL, ylim=YL, expand=FALSE)
  ggsave(file, q, width=13, height=6, dpi=220); cat("  wrote", file, "\n") }
zoom("Drosera_capensis chr15_collapsed",  "DR/fig/DR25_zoom_capensis_chr15.png")
zoom("Drosera_regia chr5_collapsed",      "DR/fig/DR25_zoom_regia_chr5.png")
zoom("Drosera_scorpioides chr1_hap1",     "DR/fig/DR25_zoom_scorpioides_chr1.png")
cat("\n  DONE.\n")
