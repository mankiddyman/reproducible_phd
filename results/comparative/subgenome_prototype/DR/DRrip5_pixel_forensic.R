#!/usr/bin/env Rscript
# ============================================================================
# DR23 — what is ACTUALLY rendered around the bar on capensis chr15_collapsed
# Gene-rank coordinates (DR20 map, exact). Dumps the pixel neighbourhood as
# numbers, then renders 4 zoom crops: current + 3 candidate fixes.
# OUT DR/fig/DR23_zoom_{current,mask,boxfill,subgcol}.png, DR/out/DR23_*.csv
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({ library(GENESPACE); library(data.table); library(ggplot2) })
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
GS <- file.path(dirname(getwd()), "genespace")
KEY <- "Drosera_capensis chr15_collapsed"
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }

hr("1. MAP (gene rank, DR20)")
e <- new.env(); load(file.path(GS,"riparian","Nepenthes_gracilis_geneOrder_rSourceData.rda"), envir=e)
P <- e$srcd$ggplotObj
BR    <- as.data.frame(P$layers[[1]]$data)
CHROM <- as.data.frame(e$srcd$sourceData$chromosomes); CHROM$key <- paste(CHROM$genome, CHROM$chr)
BLK   <- as.data.frame(e$srcd$sourceData$blocks)
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
G <- read.csv("DR/out/DR02_propagated_blocks.csv", stringsAsFactors=FALSE)
G$key <- paste(G$genome, G$chr); G <- G[G$label %in% c("A","B"), ]
m <- match(paste(BED$key,BED$id), paste(G$key,G$gene)); BED$lab <- G$label[m]; BED$reg <- G$region[m]
PAIRS <- read.csv("fractionation_by_chrpair.csv", stringsAsFactors=FALSE); s1 <- PAIRS$retained_more==PAIRS$chrA
DIOSIDE <- c(setNames(ifelse(s1,"A","B"),PAIRS$chrA), setNames(ifelse(s1,"B","A"),PAIRS$chrB))
di <- BED$genome=="Dionaea_muscipula"; BED$lab[di] <- unname(DIOSIDE[BED$chr[di]]); BED$reg[di] <- "Dionaea"
cat("  map built; genes labelled:", sum(!is.na(BED$lab)), "\n")

hr("2. TRACTS + CLASSES")
CHR <- data.frame(key=CHROM$key, genome=CHROM$genome, chr=CHROM$chr,
                  x0=CHROM$x1, x1=CHROM$x2, by1=pmin(CHROM$y1,CHROM$y2),
                  by2=pmax(CHROM$y1,CHROM$y2), y=(CHROM$y1+CHROM$y2)/2, stringsAsFactors=FALSE)
LG <- BED[!is.na(BED$lab) & !di, ]
TR <- do.call(rbind, lapply(split(LG, LG$key), function(z){ z <- z[order(z$px), ]
  r <- rle(paste(z$reg,z$lab)); en <- cumsum(r$lengths); st <- en-r$lengths+1
  do.call(rbind, lapply(seq_along(r$lengths), function(j){ w <- z[st[j]:en[j], ]
    data.frame(key=w$key[1], genome=w$genome[1], chr=w$chr[1], region=w$reg[1], label=w$lab[1],
               xa=min(w$px), xb=max(w$px), n=nrow(w), stringsAsFactors=FALSE) })) }))
TRall <- TR; TR <- TR[TR$n >= 10, ]
DC <- CHR[CHR$genome=="Dionaea_muscipula", ]
IV <- rbind(merge(TR[,c("key","genome","chr","region","label","xa","xb","n")],
                  CHR[,c("key","y","by1","by2")], by="key"),
            data.frame(key=DC$key, genome=DC$genome, chr=DC$chr, region="Dionaea",
                       label=unname(DIOSIDE[DC$chr]), xa=DC$x0, xb=DC$x1, n=NA,
                       y=DC$y, by1=DC$by1, by2=DC$by2, stringsAsFactors=FALSE))
EP <- as.data.frame(as.data.table(BR)[, .(xlo=min(x), xhi=max(x)), by=.(blkID, yrow=round(y,6))])
SD <- rbind(
  data.frame(blkID=BLK$blkID, yrow=round(as.numeric(BLK$y1),6), key=paste(BLK$genome1,BLK$chr1),
             bpa=BLK$startBp1, bpb=BLK$endBp1, g1=BLK$firstGene1, g2=BLK$lastGene1, stringsAsFactors=FALSE),
  data.frame(blkID=BLK$blkID, yrow=round(as.numeric(BLK$y2),6), key=paste(BLK$genome2,BLK$chr2),
             bpa=BLK$startBp2, bpb=BLK$endBp2, g1=BLK$firstGene2, g2=BLK$lastGene2, stringsAsFactors=FALSE))
SD <- merge(SD, EP, by=c("blkID","yrow"))
kk <- paste(BED$key, BED$ofID)
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
ovl <- function(x0,x1,y,cls){ k <- which(abs(IV$y-y)<0.25 & IV$xb>x0 & IV$xa<x1)
  if(!length(k)) return(0); w <- pmin(IV$xb[k],x1)-pmax(IV$xa[k],x0); sum(w[IV$label[k]!=cls])/max(x1-x0,1e-9) }
D <- W[W$cls %in% c("A","B"), ]
D$f1 <- mapply(ovl,D$xlo1,D$xhi1,D$yrow1,D$cls); D$f2 <- mapply(ovl,D$xlo2,D$xhi2,D$yrow2,D$cls)
KEEP <- D$blkID[pmax(D$f1,D$f2) <= 0]
BR$cls <- W$cls[match(BR$blkID,W$blkID)]; BR$cls[is.na(BR$cls) | !(BR$blkID %in% KEEP)] <- "no call"
lum <- function(h_,s_,v_){h <- rgb2hsv(col2rgb(h_)); hsv(h["h",],pmin(1,h["s",]*s_),pmin(1,h["v",]*v_))}
cols <- unique(BR$color); Ac <- setNames(lum(cols,0.55,1.10),cols); Bc <- setNames(lum(cols,1.00,0.88),cols)
BR$fillcol <- ifelse(BR$cls=="A", Ac[BR$color], Bc[BR$color])
BARA <- "#1D9E75"; BARB <- "#D85A30"
cat("  drawn braids:", length(KEEP), "| classes:\n"); print(table(BR$cls[!duplicated(BR$blkID)]))

hr("3. PALETTE COLLISION — does any braid fill read as a BAR colour?")
hu <- function(h) as.numeric(rgb2hsv(col2rgb(h))["h",])*360
near <- function(h){ d <- function(a,b){ p <- col2rgb(a); q <- col2rgb(b); sqrt(sum((p-q)^2)) }
  c(dA=d(h,BARA), dB=d(h,BARB)) }
PL <- data.frame(region_color=cols, A_fill=unname(Ac[cols]), B_fill=unname(Bc[cols]), stringsAsFactors=FALSE)
PL$hue_A <- round(hu(PL$A_fill)); PL$hue_B <- round(hu(PL$B_fill))
dA <- do.call(rbind, lapply(PL$A_fill, near)); dB <- do.call(rbind, lapply(PL$B_fill, near))
PL$A_dist_green <- dA[,1]; PL$A_dist_orange <- dA[,2]
PL$B_dist_green <- dB[,1]; PL$B_dist_orange <- dB[,2]
PL$A_reads <- ifelse(PL$A_dist_green < PL$A_dist_orange, "GREEN(=A bar)", "orange(=B bar)")
PL$B_reads <- ifelse(PL$B_dist_green < PL$B_dist_orange, "GREEN(=A bar)", "orange(=B bar)")
rg <- BR[!duplicated(BR$blkID), c("blkID","color","cls")]
rg$region <- sub("^.*_[0-9]+ ","",rg$blkID)
PL$regions <- sapply(PL$region_color, function(c0) paste(unique(rg$region[rg$color==c0]), collapse=","))
print(PL[, c("regions","region_color","A_fill","hue_A","A_reads","B_fill","hue_B","B_reads")], row.names=FALSE)
cat(sprintf("\n  B braids (vivid) whose fill is nearer the GREEN bar colour: %d of %d hues\n",
            sum(PL$B_reads=="GREEN(=A bar)"), nrow(PL)))
cat(sprintf("  A braids (muted) whose fill is nearer the GREEN bar colour: %d of %d hues\n",
            sum(PL$A_reads=="GREEN(=A bar)"), nrow(PL)))
write.csv(PL, "DR/out/DR23_palette.csv", row.names=FALSE)

hr(paste("4. PIXEL NEIGHBOURHOOD —", KEY))
ci <- which(CHR$key==KEY); Y <- CHR$y[ci]; X0 <- CHR$x0[ci]; X1 <- CHR$x1[ci]
BTOP <- Y-0.05; BBOT <- Y-0.11
cat(sprintf("  box x %.0f..%.0f  y %.4f..%.4f | bar band y %.4f..%.4f\n",
            X0, X1, CHR$by1[ci], CHR$by2[ci], BBOT, BTOP))
cat(sprintf("  braids attach at y = %.4f -- that is %s the bar band\n",
            CHR$by1[ci], ifelse(CHR$by1[ci] > BBOT & CHR$by1[ci] < BTOP, "INSIDE", "outside")))
b <- IV[IV$key==KEY, ]; b <- b[order(b$xa), ]
b$pct_start <- round(100*(b$xa-X0)/(X1-X0),1); b$pct_end <- round(100*(b$xb-X0)/(X1-X0),1)
b$colour <- ifelse(b$label=="A", BARA, BARB)
cat("\n  BARS DRAWN:\n"); print(b[,c("region","label","n","pct_start","pct_end","colour")], row.names=FALSE)
gp <- data.frame(from=c(X0,head(b$xb,-1)+1), to=c(b$xa-1,X1))
gp <- gp[gp$to > gp$from+2, ]
if (nrow(gp)) { gp$pct_from <- round(100*(gp$from-X0)/(X1-X0),1); gp$pct_to <- round(100*(gp$to-X0)/(X1-X0),1)
  gp$width_genes <- round(gp$to-gp$from)
  cat("\n  GAPS with NO bar (braids show through here):\n"); print(gp, row.names=FALSE) } else
  cat("\n  no gaps in bar coverage\n")
drop <- TRall[TRall$key==KEY & TRall$n < 10, ]
cat(sprintf("  tracts dropped for n<10 on this chromosome: %d (%s)\n", nrow(drop),
            paste(paste0(drop$region,"/",drop$label,":",drop$n), collapse=" ")))

BRd <- as.data.table(BR[BR$cls %in% c("A","B"), ])
att <- BRd[abs(y-CHR$by1[ci])<1e-9 & x>=X0 & x<=X1, .(x0=min(x), x1=max(x)), by=.(blkID,cls,fillcol)]
cat(sprintf("\n  BRAIDS ATTACHING here: %d\n", nrow(att)))
if (nrow(att)) { att$region <- sub("^.*_[0-9]+ ","",att$blkID)
  att$pct <- round(100*((att$x0+att$x1)/2-X0)/(X1-X0),1)
  att$bar_under <- sapply((att$x0+att$x1)/2, function(v){ k <- which(b$xa<=v & b$xb>=v)
    if(length(k)) paste0(b$label[k[1]],"(",b$region[k[1]],")") else "GAP" })
  att$reads_as <- ifelse(sapply(att$fillcol,function(h) near(h)["dA"] < near(h)["dB"]),"green","orange")
  print(att[order(att$pct), .(region,cls,fillcol,reads_as,pct,bar_under)], row.names=FALSE) }

cr <- BRd[y<=BTOP & y>=BBOT & x>=X0 & x<=X1, .(x0=min(x),x1=max(x)), by=.(blkID,cls,fillcol)]
cr <- cr[!(blkID %in% att$blkID)]
bel <- BRd[y<BBOT & y>=BBOT-0.08 & x>=X0 & x<=X1, .(x0=min(x),x1=max(x)), by=.(blkID,cls,fillcol)]
cat(sprintf("\n  FOREIGN braids crossing INSIDE the bar band: %d (hidden -- bars drawn on top)\n", nrow(cr)))
cat(sprintf("  braids visible in the 0.08 strip JUST BELOW the bar: %d\n", nrow(bel)))
if (nrow(bel)) { bel$reads_as <- ifelse(sapply(bel$fillcol,function(h) near(h)["dA"] < near(h)["dB"]),"green","orange")
  bel$pct <- round(100*((bel$x0+bel$x1)/2-X0)/(X1-X0),1)
  bel$bar_above <- sapply((bel$x0+bel$x1)/2, function(v){ k <- which(b$xa<=v & b$xb>=v)
    if(length(k)) b$label[k[1]] else "GAP" })
  bel$CLASH <- (bel$reads_as=="green" & bel$bar_above=="B") | (bel$reads_as=="orange" & bel$bar_above=="A")
  print(bel[order(bel$pct), .(cls,fillcol,reads_as,pct,bar_above,CLASH)], row.names=FALSE)
  cat(sprintf("\n  >>> GREEN-LOOKING pixels sitting directly below an ORANGE bar: %d\n",
              sum(bel$reads_as=="green" & bel$bar_above=="B")))
  write.csv(bel, "DR/out/DR23_below_bar.csv", row.names=FALSE) }

hr("5. FOUR ZOOM CROPS")
XL <- c(X0-40, X1+40); YL <- c(Y-0.75, Y+0.14)
keepb <- unique(BR$blkID[BR$x>=XL[1] & BR$x<=XL[2] & BR$y>=YL[1] & BR$y<=YL[2]])
BZ <- BR[BR$blkID %in% keepb, ]
mkl <- function(d,al,fc="fillcol") layer(geom="polygon", stat="identity", position="identity", data=d,
  mapping=aes(x=x,y=y,group=blkID,fill=.data[[fc]]), params=list(alpha=al,colour=NA,linewidth=0), show.legend=FALSE)
rects <- function(d,fillvar) layer(geom="rect", stat="identity", position="identity", data=d,
  mapping=aes(xmin=xa,xmax=xb,ymin=ybot,ymax=ytop,fill=.data[[fillvar]]), params=list(colour=NA), show.legend=FALSE)
render <- function(braid, extras, pal, file, sub){
  q <- P; q$layers <- c(braid, P$layers[-1], extras)
  q$scales$scales <- q$scales$scales[!vapply(q$scales$scales,function(s)"fill"%in%s$aesthetics,logical(1))]
  q <- q + scale_fill_manual(values=pal, guide="none") +
    coord_cartesian(xlim=XL, ylim=YL, expand=FALSE) + theme_minimal(11) +
    theme(panel.background=element_rect(fill="grey96",colour=NA),
          plot.background=element_rect(fill="grey96",colour=NA), panel.grid=element_blank(),
          axis.text=element_blank(), axis.ticks=element_blank(),
          plot.subtitle=element_text(size=9)) +
    labs(title=KEY, subtitle=sub, x=NULL, y=NULL)
  ggsave(file, q, width=13, height=5, dpi=220); cat("  wrote", file, "\n") }

IVz <- IV; IVz$barcol <- ifelse(IVz$label=="A",BARA,BARB)
IVz$ytop <- IVz$y-0.05; IVz$ybot <- IVz$y-0.11
PALb <- c(setNames(unique(BZ$fillcol),unique(BZ$fillcol)), setNames(c(BARA,BARB),c(BARA,BARB)))
bl <- list(mkl(BZ[BZ$cls=="A",],0.55), mkl(BZ[BZ$cls=="B",],0.88))
render(bl, list(rects(IVz,"barcol")), PALb, "DR/fig/DR23_zoom_current.png",
       "CURRENT: bar at y-0.05..-0.11; braids attach inside it and emerge below")

MK <- unique(IVz[,c("key","y")]); MK <- merge(MK, CHR[,c("key","x0","x1")], by="key")
MK$xa <- MK$x0-20; MK$xb <- MK$x1+20; MK$ytop <- MK$y-0.11; MK$ybot <- MK$y-0.20; MK$barcol <- "grey96"
PALm <- c(PALb, c(grey96="grey96"))
render(bl, list(rects(MK,"barcol"), rects(IVz,"barcol")), PALm, "DR/fig/DR23_zoom_mask.png",
       "FIX 1: background strip below each bar, so nothing foreign touches it")

h <- IVz$by2-IVz$by1; IVf <- IVz; IVf$ybot <- IVz$by1+0.18*h; IVf$ytop <- IVz$by2-0.18*h
render(bl, list(rects(IVf,"barcol")), PALb, "DR/fig/DR23_zoom_boxfill.png",
       "FIX 2: tracts filled INTO the chromosome box -- no braid can be adjacent")

BZ2 <- BZ; BZ2$sgcol <- ifelse(BZ2$cls=="A",BARA,BARB)
PALs <- c(setNames(c(BARA,BARB),c(BARA,BARB)))
render(list(mkl(BZ2[BZ2$cls=="A",],0.75,"sgcol"), mkl(BZ2[BZ2$cls=="B",],0.90,"sgcol")),
       list(rects(IVz,"barcol")), PALs, "DR/fig/DR23_zoom_subgcol.png",
       "FIX 3: braids coloured by SUBGENOME, same green/orange as the bars (region hue dropped)")
