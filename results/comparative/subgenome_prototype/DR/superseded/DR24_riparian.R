#!/usr/bin/env Rscript
# ============================================================================
# DR24 — PRODUCTION riparian: gene rank, tracts filled INTO the chromosome box
#
# Coordinate map (DR18/DR19/DR20, proven to 0.000 plot units over 2556 sides):
#     plot_ord = rank(ord) within chromosome - offset(chromosome)
#     x        = chromosome$x1 + plot_ord           [no chromosome is flipped]
#
# The defect this fixes: bars were drawn at y-0.05..y-0.11 while braids attach
# at the box edge (capensis: 6.9385, INSIDE that band). The attachment was
# painted over by the bar; the first visible braid pixels were 0.048 lower,
# after sideways drift, where 20% of ribbons pass over a wrong-coloured bar.
# Here tracts are strips at the TOP and BOTTOM EDGES of the box -- exactly
# where braids attach -- so a braid touches its own subgenome colour or
# nothing. Middle of the box stays white for the chromosome label.
#
# OUT DR/fig/DR24_riparian_AB.{pdf,png}
#     DR/fig/DR24_zoom_capensis_chr15.png, DR/fig/DR24_zoom_regia_chr5.png
#     DR/out/DR24_dropped.csv, DR/out/DR24_tracts.csv
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({ library(GENESPACE); library(data.table); library(ggplot2) })
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
GS <- file.path(dirname(getwd()), "genespace")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }
BARA <- "#1D9E75"; BARB <- "#D85A30"

hr("1. COORDINATE MAP")
e <- new.env(); load(file.path(GS,"riparian","Nepenthes_gracilis_geneOrder_rSourceData.rda"), envir=e)
P <- e$srcd$ggplotObj
BR    <- as.data.frame(P$layers[[1]]$data)
CHROM <- as.data.frame(e$srcd$sourceData$chromosomes); CHROM$key <- paste(CHROM$genome, CHROM$chr)
BLK   <- as.data.frame(e$srcd$sourceData$blocks)
cat(sprintf("  saved plot layers: %d | braid rows %d | blocks %d\n",
            length(P$layers), nrow(BR), nrow(BLK)))
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
write.csv(IV, "DR/out/DR24_tracts.csv", row.names=FALSE)

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
                          "DR/out/DR24_dropped.csv", row.names=FALSE)
BR$cls <- W$cls[match(BR$blkID,W$blkID)]; BR$cls[is.na(BR$cls) | !(BR$blkID %in% KEEP)] <- "no call"
chk <- D[D$blkID %in% KEEP, ]
if (nrow(chk) && max(chk$fmax) > 0) die("enforcement failed")
cat("  GUARANTEE: every drawn attachment lies wholly on strips of its own subgenome.\n")

hr("6. DRAW — tracts as strips at the box EDGES")
lum <- function(h_,s_,v_){h <- rgb2hsv(col2rgb(h_)); hsv(h["h",],pmin(1,h["s",]*s_),pmin(1,h["v",]*v_))}
cols <- unique(BR$color); Ac <- setNames(lum(cols,0.55,1.10),cols); Bc <- setNames(lum(cols,1.00,0.88),cols)
BR$fillcol <- ifelse(BR$cls=="A", Ac[BR$color], Bc[BR$color])
h <- IV$by2 - IV$by1; FR <- 0.30                      # strip height, fraction of box
ST <- rbind(
  transform(IV, ymin = by2 - FR*h, ymax = by2),       # top edge
  transform(IV, ymin = by1,        ymax = by1 + FR*h))# bottom edge
ST$barcol <- ifelse(ST$label=="A", BARA, BARB)
SEP <- rbind(transform(IV, ymin=by2-0.012, ymax=by2, xa=x0, xb=x1),
             transform(IV, ymin=by1, ymax=by1+0.012, xa=x0, xb=x1))
SEP <- SEP[!duplicated(paste(SEP$key, SEP$ymin)), ]; SEP$barcol <- "#FFFFFF"
PAL <- c(setNames(unique(BR$fillcol),unique(BR$fillcol)),
         setNames(c(BARA,BARB,"#FFFFFF"), c(BARA,BARB,"#FFFFFF")))
mkl <- function(d,al) layer(geom="polygon", stat="identity", position="identity", data=d,
  mapping=aes(x=x,y=y,group=blkID,fill=fillcol), params=list(alpha=al,colour=NA,linewidth=0), show.legend=FALSE)
rects <- function(d) layer(geom="rect", stat="identity", position="identity", data=d,
  mapping=aes(xmin=xa,xmax=xb,ymin=ymin,ymax=ymax,fill=barcol), params=list(colour=NA), show.legend=FALSE)
nL <- length(P$layers)
build <- function(braid) { q <- P
  q$layers <- c(braid, P$layers[2], list(rects(ST), rects(SEP)),
                if (nL >= 3) P$layers[3:nL] else NULL)
  q$scales$scales <- q$scales$scales[!vapply(q$scales$scales,function(s)"fill"%in%s$aesthetics,logical(1))]
  q + scale_fill_manual(values=PAL, guide="none") }
sty <- function(q, ttl, sub) q + theme_minimal(11) + theme(
  panel.background=element_rect(fill="grey96",colour=NA), plot.background=element_rect(fill="grey96",colour=NA),
  panel.grid=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(),
  axis.text.y=element_text(face="italic",size=11), plot.title=element_text(face="bold",size=14),
  plot.subtitle=element_text(size=9,colour="grey30")) +
  labs(title=ttl, subtitle=sub, x="Chromosomes scaled by GENE RANK -- physical spacing is not preserved", y=NULL)
BL <- list(mkl(BR[BR$cls=="A",],0.55), mkl(BR[BR$cls=="B",],0.88))
SUB <- paste0("Coloured strips along each chromosome edge are the phased subgenome tracts: green = A, orange = B. ",
  "A braid is drawn only where\nboth ends are >=95% one subgenome among the block's own genes AND attach wholly ",
  "within strips of that subgenome.\nMUTED = A, VIVID = B; hue = ancestral region. Braid hue is independent of the ",
  "green/orange subgenome coding.")
p <- sty(build(BL), "Droseraceae riparian, phased by subgenome", SUB)
ggsave("DR/fig/DR24_riparian_AB.pdf", p, width=16, height=9)
ggsave("DR/fig/DR24_riparian_AB.png", p, width=16, height=9, dpi=200)
cat(sprintf("  wrote DR/fig/DR24_riparian_AB.{pdf,png} | %d braids drawn\n", length(KEEP)))

hr("7. ZOOM VERIFICATION")
zoom <- function(key, file){ i <- which(CHR$key==key); if(!length(i)) return(invisible())
  XL <- c(CHR$x0[i]-60, CHR$x1[i]+60); YL <- c(CHR$y[i]-0.62, CHR$y[i]+0.62)
  kb <- unique(BR$blkID[BR$x>=XL[1] & BR$x<=XL[2] & abs(BR$y-CHR$y[i])<0.62])
  BZ <- BR[BR$blkID %in% kb, ]
  q <- build(list(mkl(BZ[BZ$cls=="A",],0.55), mkl(BZ[BZ$cls=="B",],0.88)))
  q <- sty(q, key, "zoom -- every braid should touch a strip of its own colour") +
       coord_cartesian(xlim=XL, ylim=YL, expand=FALSE)
  ggsave(file, q, width=13, height=6, dpi=220); cat("  wrote", file, "\n")
  b <- IV[IV$key==key, ]; b <- b[order(b$xa), ]
  b$pct_start <- round(100*(b$xa-CHR$x0[i])/(CHR$x1[i]-CHR$x0[i]),1)
  b$pct_end   <- round(100*(b$xb-CHR$x0[i])/(CHR$x1[i]-CHR$x0[i]),1)
  print(b[,c("region","label","n","pct_start","pct_end")], row.names=FALSE) }
zoom("Drosera_capensis chr15_collapsed", "DR/fig/DR24_zoom_capensis_chr15.png")
zoom("Drosera_regia chr5_collapsed",     "DR/fig/DR24_zoom_regia_chr5.png")
cat("\n  DONE.\n")
