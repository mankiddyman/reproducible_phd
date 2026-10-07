#!/usr/bin/env Rscript
# ============================================================================
# DR20 — DR19's proven coordinate map + an audit that measures what the EYE sees
#
# DR19 fixed geometry (gene->x exact, 0.000 plot units over 2556 block sides)
# but inherited DR17's two labelling defects:
#   endlab()    scored a braid end using only genes of the BLOCK's region
#   barlab_at() ignored any bar whose region differed, and tested the MIDPOINT
# Both make a vivid B braid over a green A bar invisible to the audit.
#
# Here:
#   - a braid end is labelled from the BLOCK'S OWN GENES (rank interval between
#     firstGene and lastGene), not from a window+region proxy
#   - the audit is region-agnostic and span-weighted
#   - any end overlapping a wrong-coloured bar is DROPPED before drawing, so
#     the defect cannot appear in the figure by construction
#
# OUT  DR/fig/DR20_riparian_AB.{pdf,png}, DR/out/DR20_dropped.csv,
#      DR/out/DR20_overlap_audit.csv
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({ library(GENESPACE); library(data.table); library(ggplot2) })
setwd(Sys.getenv("SUBG_BASE", getwd()))
GS <- file.path(dirname(getwd()), "genespace")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }

hr("1. LOAD + CANONICAL GENE -> PLOT X  (DR19 section 2, unchanged)")
e <- new.env()
load(file.path(GS,"riparian","Nepenthes_gracilis_geneOrder_rSourceData.rda"), envir=e)
P <- e$srcd$ggplotObj
BR    <- as.data.frame(P$layers[[1]]$data)
CHROM <- as.data.frame(e$srcd$sourceData$chromosomes); CHROM$key <- paste(CHROM$genome, CHROM$chr)
BLK   <- as.data.frame(e$srcd$sourceData$blocks)
BED <- as.data.frame(fread(file.path(GS,"results","combBed.txt")))
BED$key <- paste(BED$genome, BED$chr); BED <- BED[BED$key %in% CHROM$key, ]
BED <- BED[order(BED$key, BED$ord), ]
BED$r <- ave(BED$ord, BED$key, FUN=function(z) z - min(z) + 1L)
AN <- unique(rbind(
  data.frame(key=paste(BLK$genome1,BLK$chr1), gene=BLK$firstGene1, pord=BLK$startOrd1),
  data.frame(key=paste(BLK$genome1,BLK$chr1), gene=BLK$lastGene1,  pord=BLK$endOrd1),
  data.frame(key=paste(BLK$genome2,BLK$chr2), gene=BLK$firstGene2, pord=BLK$startOrd2),
  data.frame(key=paste(BLK$genome2,BLK$chr2), gene=BLK$lastGene2,  pord=BLK$endOrd2),
  stringsAsFactors=FALSE))
AN$r <- BED$r[match(paste(AN$key,AN$gene), paste(BED$key,BED$ofID))]; AN <- AN[!is.na(AN$r), ]
AN$d <- AN$r - AN$pord
if (max(tapply(AN$d, AN$key, function(z) diff(range(z)))) != 0) die("offset not constant")
OFF <- tapply(AN$d, AN$key, function(z) z[1])
BED$cx1 <- CHROM$x1[match(BED$key,CHROM$key)]; BED$clen <- CHROM$length[match(BED$key,CHROM$key)]
BED$px  <- BED$cx1 + pmin(pmax(BED$r - OFF[BED$key], 1), BED$clen)
EP <- as.data.frame(as.data.table(BR)[, .(xlo=min(x), xhi=max(x)), by=.(blkID, yrow=round(y,6))])
kk <- paste(BED$key, BED$ofID)
SD <- rbind(
  data.frame(blkID=BLK$blkID, yrow=round(as.numeric(BLK$y1),6), key=paste(BLK$genome1,BLK$chr1),
             ga=BLK$firstGene1, gb=BLK$lastGene1, stringsAsFactors=FALSE),
  data.frame(blkID=BLK$blkID, yrow=round(as.numeric(BLK$y2),6), key=paste(BLK$genome2,BLK$chr2),
             ga=BLK$firstGene2, gb=BLK$lastGene2, stringsAsFactors=FALSE))
SD <- merge(SD, EP, by=c("blkID","yrow"))
SD$xa <- BED$px[match(paste(SD$key,SD$ga), kk)]; SD$xb <- BED$px[match(paste(SD$key,SD$gb), kk)]
SD <- SD[!is.na(SD$xa) & !is.na(SD$xb), ]
err <- max(pmax(abs(pmin(SD$xa,SD$xb)-SD$xlo), abs(pmax(SD$xa,SD$xb)-SD$xhi)))
cat(sprintf("  block sides %d | max coordinate error %.4g plot units\n", nrow(SD), err))
if (err > 1e-6) die("coordinate map regressed")

hr("2. GENE LABELS  (A/B on every gene, one table)")
G <- read.csv("DR/out/DR02_propagated_blocks.csv", stringsAsFactors=FALSE)
gcol <- intersect(c("gene","id","ofID"), names(G))[1]
G$key <- paste(G$genome, G$chr); G <- G[!is.na(G$label) & G$label %in% c("A","B"), ]
m1 <- match(paste(BED$key,BED$id), paste(G$key,G[[gcol]]))
m2 <- match(kk,                    paste(G$key,G[[gcol]]))
m  <- if (sum(!is.na(m1)) >= sum(!is.na(m2))) m1 else m2
BED$lab <- G$label[m]; BED$reg <- G$region[m]
PAIRS <- read.csv("fractionation_by_chrpair.csv", stringsAsFactors=FALSE)
s1 <- PAIRS$retained_more == PAIRS$chrA
DIOSIDE <- c(setNames(ifelse(s1,"A","B"), PAIRS$chrA), setNames(ifelse(s1,"B","A"), PAIRS$chrB))
di <- BED$genome == "Dionaea_muscipula"
BED$lab[di] <- unname(DIOSIDE[BED$chr[di]]); BED$reg[di] <- "Dionaea"
cat(sprintf("  labelled genes: %d of %d | Drosera %d | Dionaea %d | Nepenthes %d (outgroup, unlabelled)\n",
    sum(!is.na(BED$lab)), nrow(BED), sum(!is.na(BED$lab) & !di),
    sum(!is.na(BED$lab) & di), sum(BED$genome=="Nepenthes_gracilis")))

hr("3. TRACTS = THE BARS")
CHR <- data.frame(genome=CHROM$genome, chr=CHROM$chr, x0=CHROM$x1, x1=CHROM$x2,
                  y=(CHROM$y1+CHROM$y2)/2, key=CHROM$key, stringsAsFactors=FALSE)
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
TR$y <- CHR$y[match(paste(TR$genome,TR$chr), paste(CHR$genome,CHR$chr))]
DC <- CHR[CHR$genome=="Dionaea_muscipula", ]; DC$label <- unname(DIOSIDE[DC$chr])
IV <- rbind(TR[, c("genome","chr","region","label","xa","xb","y")],
            data.frame(genome=DC$genome, chr=DC$chr, region="Dionaea", label=DC$label,
                       xa=DC$x0, xb=DC$x1, y=DC$y, stringsAsFactors=FALSE))
IV$ytop <- IV$y - 0.05; IV$ybot <- IV$y - 0.11
cat(sprintf("  tracts %d (%d A, %d B) | chromosomes with both labels %d | bars drawn %d\n",
    nrow(TR), sum(TR$label=="A"), sum(TR$label=="B"), .both, nrow(IV)))
ck <- do.call(rbind, lapply(split(IV, paste(IV$genome,IV$chr)), function(z){
  z <- z[order(z$xa), ]; data.frame(ov=sum(head(z$xb,-1) > tail(z$xa,-1))) }))
cat(sprintf("  overlapping bar pairs on the same chromosome: %d (must be 0)\n", sum(ck$ov)))
if (sum(ck$ov) > 0) die("bars overlap -- a point would have two colours")

hr("4. BRAID ENDS LABELLED FROM THE BLOCK'S OWN GENES")
DT <- as.data.table(BED[, c("key","r","px","lab")]); setkey(DT, key, r)
SD$ra <- BED$r[match(paste(SD$key,SD$ga), kk)]; SD$rb <- BED$r[match(paste(SD$key,SD$gb), kk)]
SD$lo <- pmin(SD$ra, SD$rb); SD$hi <- pmax(SD$ra, SD$rb)
cmp <- t(vapply(seq_len(nrow(SD)), function(i){
  z <- DT[.(SD$key[i])][r >= SD$lo[i] & r <= SD$hi[i], lab]
  z <- z[!is.na(z)]; if (!length(z)) return(c(0,0))
  c(sum(z=="A"), sum(z=="B")) }, numeric(2)))
SD$nA <- cmp[,1]; SD$nB <- cmp[,2]; SD$nTot <- SD$nA + SD$nB
SD$pur <- ifelse(SD$nTot > 0, pmax(SD$nA,SD$nB)/SD$nTot, NA_real_)
SD$L   <- ifelse(is.na(SD$pur) | SD$pur < 0.95, NA_character_,
                 ifelse(SD$nA >= SD$nB, "A", "B"))
cat(sprintf("  block sides with >=1 labelled gene of their own: %d of %d\n",
    sum(SD$nTot > 0), nrow(SD)))
cat("  genes per block side:\n"); print(summary(SD$nTot))
cat("  purity of the block's own genes:\n"); print(summary(SD$pur))

W <- merge(SD[SD$yrow==round(as.numeric(BLK$y1),6)[match(SD$blkID,BLK$blkID)],
              c("blkID","key","yrow","xlo","xhi","L","pur","nTot")],
           SD[, c("blkID","yrow","key","xlo","xhi","L","pur","nTot")],
           by="blkID", suffixes=c("1","2"))
W <- W[W$yrow1 != W$yrow2, ]; W <- W[!duplicated(W$blkID), ]
nep <- function(k) startsWith(k, "Nepenthes_gracilis ")
W$cls <- ifelse(!is.na(W$L1) & !is.na(W$L2), ifelse(W$L1==W$L2, W$L1, "discordant"),
          ifelse(is.na(W$L1) & nep(W$key1) & !is.na(W$L2), W$L2,
          ifelse(is.na(W$L2) & nep(W$key2) & !is.na(W$L1), W$L1, "no call")))
cat("\n  braids by class (from the block's own genes):\n"); print(table(W$cls))

hr("5. HONEST AUDIT — region-agnostic, span-weighted overlap with the bars")
ovl <- function(x0, x1, y){
  k <- which(abs(IV$y - y) < 0.25 & IV$xb > x0 & IV$xa < x1)
  if (!length(k)) return(c(0,0))
  w <- pmin(IV$xb[k],x1) - pmax(IV$xa[k],x0)
  c(sum(w[IV$label[k]=="A"]), sum(w[IV$label[k]=="B"])) }
DR <- W[W$cls %in% c("A","B"), ]
o1 <- t(mapply(ovl, DR$xlo1, DR$xhi1, DR$yrow1))
o2 <- t(mapply(ovl, DR$xlo2, DR$xhi2, DR$yrow2))
DR$span1 <- DR$xhi1-DR$xlo1; DR$span2 <- DR$xhi2-DR$xlo2
DR$wrong1 <- ifelse(DR$cls=="A", o1[,2], o1[,1])
DR$wrong2 <- ifelse(DR$cls=="A", o2[,2], o2[,1])
DR$f1 <- DR$wrong1/pmax(DR$span1,1); DR$f2 <- DR$wrong2/pmax(DR$span2,1)
DR$fmax <- pmax(DR$f1, DR$f2); DR$wmax <- pmax(DR$wrong1, DR$wrong2)
cat(sprintf("  candidate braids: %d\n", nrow(DR)))
cat("  fraction of an end's span sitting on a WRONG-coloured bar:\n")
print(summary(DR$fmax))
for (th in c(0, 0.01, 0.05, 0.10, 0.25)) cat(sprintf(
  "    ends exceeding %5.0f%%: %4d braids (%5.1f%%) -- would remain: %d\n",
  100*th, sum(DR$fmax > th), 100*mean(DR$fmax > th), sum(DR$fmax <= th)))
cat("\n  worst 12 offenders:\n")
print(head(DR[order(-DR$fmax), c("key1","key2","cls","fmax","wmax","pur1","pur2","nTot1","nTot2")], 12))
cat("\n  --- Drosera_capensis chr15_collapsed, every candidate braid ---\n")
cc <- DR[DR$key1=="Drosera_capensis chr15_collapsed" | DR$key2=="Drosera_capensis chr15_collapsed", ]
print(cc[order(-cc$fmax), c("blkID","cls","f1","f2","wrong1","wrong2","span1","span2","pur1","pur2")])
write.csv(DR[, c("blkID","key1","key2","cls","f1","f2","wrong1","wrong2","span1","span2",
                 "pur1","pur2","nTot1","nTot2")], "DR/out/DR20_overlap_audit.csv", row.names=FALSE)

hr("6. ENFORCE — drop every end that touches a wrong-coloured bar")
TOL <- 0
KEEPB <- DR$blkID[DR$fmax <= TOL]
DROP  <- DR[DR$fmax > TOL, ]
cat(sprintf("  kept %d | dropped %d (%.1f%%)\n", length(KEEPB), nrow(DROP), 100*mean(DR$fmax>TOL)))
if (nrow(DROP)) { write.csv(DROP[, c("blkID","key1","key2","cls","fmax","wmax")],
  "DR/out/DR20_dropped.csv", row.names=FALSE); cat("  wrote DR/out/DR20_dropped.csv\n") }
BR$cls <- W$cls[match(BR$blkID, W$blkID)]
BR$cls[is.na(BR$cls) | !(BR$blkID %in% KEEPB)] <- "no call"
cat("  braids actually drawn:\n"); print(table(BR$cls[!duplicated(BR$blkID)]))
res <- DR[DR$blkID %in% KEEPB, ]
if (nrow(res) && max(res$fmax) > 0) die("enforcement failed: %d ends still overlap", sum(res$fmax>0))
cat("\n  GUARANTEE: every drawn end sits entirely on bars of its own colour.\n")

hr("7. DRAW")
lum <- function(h_,s_,v_){h <- rgb2hsv(col2rgb(h_)); hsv(h["h",],pmin(1,h["s",]*s_),pmin(1,h["v",]*v_))}
cols <- unique(BR$color)
Ac <- setNames(lum(cols,0.55,1.10), cols); Bc <- setNames(lum(cols,1.00,0.88), cols)
BR$fillcol <- ifelse(BR$cls=="A", Ac[BR$color], Bc[BR$color])
IV$barcol  <- ifelse(IV$label=="A", "#1D9E75", "#D85A30")
PAL <- c(setNames(unique(BR$fillcol),unique(BR$fillcol)), setNames(unique(IV$barcol), unique(IV$barcol)))
mkl <- function(d,al) layer(geom="polygon", stat="identity", position="identity", data=d,
  mapping=aes(x=x,y=y,group=blkID,fill=fillcol),
  params=list(alpha=al,colour=NA,linewidth=0), show.legend=FALSE)
p <- P
p$layers <- c(list(mkl(BR[BR$cls=="A",],0.55), mkl(BR[BR$cls=="B",],0.88)), P$layers[-1],
  list(layer(geom="rect", stat="identity", position="identity", data=IV,
    mapping=aes(xmin=xa,xmax=xb,ymin=ybot,ymax=ytop,fill=barcol),
    params=list(colour=NA), show.legend=FALSE)))
p$scales$scales <- p$scales$scales[!vapply(p$scales$scales, function(sc) "fill" %in% sc$aesthetics, logical(1))]
p <- p + scale_fill_manual(values=PAL, guide="none") + theme_minimal(11) +
  theme(panel.background=element_rect(fill="grey96",colour=NA),
        plot.background=element_rect(fill="grey96",colour=NA), panel.grid=element_blank(),
        axis.text.x=element_blank(), axis.ticks=element_blank(),
        axis.text.y=element_text(face="italic",size=11),
        plot.title=element_text(face="bold",size=14),
        plot.subtitle=element_text(size=9,colour="grey30")) +
  labs(title="Droseraceae riparian, phased by subgenome",
       subtitle=paste0("Bars beneath each chromosome are the phased tracts: green = A, orange = B. ",
         "A braid is drawn only where both\nends are >=95% one subgenome among the block's own genes ",
         "AND lie entirely on bars of that colour. MUTED = A, VIVID = B."),
       x="Chromosomes scaled by gene rank order", y=NULL)
ggsave("DR/fig/DR20_riparian_AB.pdf", p, width=16, height=9)
ggsave("DR/fig/DR20_riparian_AB.png", p, width=16, height=9, dpi=200)
cat(sprintf("  wrote DR/fig/DR20_riparian_AB.{pdf,png} | %d braids drawn, %d dropped\n",
    length(KEEPB), nrow(DROP)))
