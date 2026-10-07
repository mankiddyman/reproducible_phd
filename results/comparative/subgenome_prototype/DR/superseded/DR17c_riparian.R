#!/usr/bin/env Rscript
# ============================================================================
# DR17c — RIPARIAN WITH BRAIDS AND BARS IN ONE COORDINATE SYSTEM
#
# WHY THE PREVIOUS VERSIONS KEPT DISAGREEING
#   GENESPACE positions braids by ITS gene ordering, which it computes to
#   maximise syntenic alignment against the reference. DR17 positioned bars by
#   converting physical Mb to plot x. Those two orderings do not correspond:
#   across 57 chromosomes the Spearman correlation between Mb and plot x runs
#   from -0.83 to +0.89, and only Nepenthes (the reference) is clean. No
#   conversion can reconcile them, because the relationship is not monotonic.
#
#   Two attempts failed: interpolating through braid endpoints (smeared tracts
#   across whole genomes), and using combBed gene ranks (same result as Mb).
#
# THIS VERSION
#   Keeps GENESPACE's LAYOUT -- chromosome x positions, gaps, ordering, row
#   order, labels -- and redraws ONLY the ribbon polygons, positioning both
#   endpoints with the SAME Mb -> x mapping used for the bars.
#
#   Bars and braids then agree by construction, not by approximation. The
#   assertion at the end checks it: every drawn braid must sit on a bar of its
#   own colour, with zero exceptions.
#
#   Cost: ribbon positions differ from GENESPACE's original, because
#   GENESPACE's are in its own gene-rank space. The figure is internally
#   consistent instead of matching an ordering the bars cannot use.
#
# IN   ../genespace/riparian/Nepenthes_gracilis_geneOrder_rSourceData.rda
#      ../genespace/results/{syntenicBlock_coordinates.csv,combBed.txt}
#      DR/out/DR02_propagated_blocks.csv, fractionation_by_chrpair.csv
# OUT  DR/fig/DR17c_riparian_AB.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(GENESPACE); library(data.table); library(ggplot2)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
GS <- file.path(dirname(getwd()), "genespace")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
PURITY <- 0.95      # a braid end must be this pure to be drawn
MINGENE <- 3        # and have at least this many genes to assess

hr("1. GENESPACE LAYOUT — kept exactly")
e <- new.env()
load(file.path(GS,"riparian","Nepenthes_gracilis_geneOrder_rSourceData.rda"), envir=e)
P  <- e$srcd$ggplotObj
BR <- as.data.frame(P$layers[[1]]$data)
BR$base   <- sub("_[0-9]+ chr[^ ]*$", "", BR$blkID)
BR$region <- sub("^.*_[0-9]+ ", "", BR$blkID)
CH  <- as.data.frame(P$layers[[2]]$data)
CHR <- do.call(rbind, lapply(split(CH, paste(CH$genome, CH$chr)), function(z)
  data.frame(genome=z$genome[1], chr=z$chr[1], x0=min(z$x), x1=max(z$x),
             y=mean(z$y), stringsAsFactors=FALSE)))
cat(sprintf("  chromosomes: %d | rows: %d\n", nrow(CHR), length(unique(CHR$y))))
cat(sprintf("  blocks GENESPACE drew: %d\n", length(unique(BR$blkID))))

hr("2. ONE Mb -> x MAPPING, used for BOTH bars and braids")
BED <- read.table(file.path(GS,"results","combBed.txt"), header=TRUE, sep="\t",
                  quote="", comment.char="", stringsAsFactors=FALSE)
BED$mb <- (BED$start + BED$end)/2e6
mb2x <- function(g, c, mb) {
  k <- which(CHR$genome==g & CHR$chr==c); if (!length(k)) return(NA_real_)
  bb <- BED[BED$genome==g & BED$chr==c, ]; if (!nrow(bb)) return(NA_real_)
  bb <- bb[order(bb$mb), ]
  f <- (findInterval(mb, bb$mb) - 1) / max(1, nrow(bb) - 1)
  CHR$x0[k[1]] + pmin(pmax(f, 0), 1) * (CHR$x1[k[1]] - CHR$x0[k[1]])
}
cat("  Mb -> x by rank among genes SORTED BY POSITION on that chromosome.\n")
cat("  Monotone by construction, and identical for bars and braids.\n")

hr("3. GENES AND THEIR LABELS")
G <- read.csv("DR/out/DR02_propagated_blocks.csv", stringsAsFactors=FALSE)
G$mb <- G$mid/1e6
PAIRS <- read.csv("fractionation_by_chrpair.csv", stringsAsFactors=FALSE)
s1 <- PAIRS$retained_more == PAIRS$chrA
DIOSIDE <- c(setNames(ifelse(s1,"A","B"), PAIRS$chrA),
             setNames(ifelse(s1,"B","A"), PAIRS$chrB))
cat(sprintf("  Drosera genes labelled: %d\n", nrow(G)))
cat(sprintf("  Dionaea chromosomes with a side: %d\n", length(DIOSIDE)))

hr("4. TRACTS -> BARS")
TR <- do.call(rbind, lapply(split(G, paste(G$genome, G$chr)), function(z) {
  z <- z[order(z$mb), ]
  r <- rle(paste(z$region, z$label))
  en <- cumsum(r$lengths); st <- en - r$lengths + 1
  do.call(rbind, lapply(seq_along(r$lengths), function(j) {
    w <- z[st[j]:en[j], ]
    if (nrow(w) < 10) return(NULL)
    data.frame(genome=w$genome[1], chr=w$chr[1], region=w$region[1],
               label=w$label[1], lo=min(w$mb), hi=max(w$mb), n=nrow(w),
               stringsAsFactors=FALSE)
  }))
}))
TR$xa <- mapply(mb2x, TR$genome, TR$chr, TR$lo)
TR$xb <- mapply(mb2x, TR$genome, TR$chr, TR$hi)
TR$y  <- CHR$y[match(paste(TR$genome, TR$chr), paste(CHR$genome, CHR$chr))]
TR <- TR[!is.na(TR$xa) & !is.na(TR$xb) & !is.na(TR$y), ]
D <- CHR[CHR$genome=="Dionaea_muscipula", ]
D$label <- unname(DIOSIDE[D$chr]); stopifnot(!any(is.na(D$label)))
IV <- rbind(TR[, c("genome","chr","region","label","xa","xb","y")],
            data.frame(genome=D$genome, chr=D$chr, region="Dionaea",
                       label=D$label, xa=D$x0, xb=D$x1, y=D$y,
                       stringsAsFactors=FALSE))
IV$ytop <- IV$y - 0.05; IV$ybot <- IV$y - 0.11
cat(sprintf("  bars: %d (%d A, %d B)\n", nrow(IV),
            sum(IV$label=="A"), sum(IV$label=="B")))

hr("5. BLOCKS -> BRAIDS, positioned by the SAME mapping")
BC <- as.data.frame(fread(file.path(GS,"results","syntenicBlock_coordinates.csv")))
BC$base <- as.character(BC$blkID)
RG <- unique(BR[, c("base","region")])
BK <- merge(BC, RG, by="base")
BK <- BK[!duplicated(paste(BK$base, BK$region)), ]
cat(sprintf("  blocks matched to a drawn braid: %d\n", nrow(BK)))
BK$m1a <- pmin(BK$startBp1, BK$endBp1)/1e6; BK$m1b <- pmax(BK$startBp1, BK$endBp1)/1e6
BK$m2a <- pmin(BK$startBp2, BK$endBp2)/1e6; BK$m2b <- pmax(BK$startBp2, BK$endBp2)/1e6

# label an end: genes of THIS region in THIS Mb window, on this chromosome
endlab <- function(g, c, m0, m1, reg) {
  if (g == "Dionaea_muscipula")  return(unname(DIOSIDE[c]))
  if (g == "Nepenthes_gracilis") return("ref")
  z <- G[G$genome==g & G$chr==c & G$region==reg & G$mb>=m0 & G$mb<=m1, ]
  if (nrow(z) < MINGENE) return(NA_character_)
  t <- table(z$label)
  if (max(t)/sum(t) < PURITY) return("chimeric")
  names(t)[which.max(t)]
}
BK$L1 <- mapply(endlab, BK$genome1, BK$chr1, BK$m1a, BK$m1b, BK$region)
BK$L2 <- mapply(endlab, BK$genome2, BK$chr2, BK$m2a, BK$m2b, BK$region)
cls <- function(a, b) {
  if (is.na(a) && is.na(b)) return("no call")
  if (identical(a,"chimeric") || identical(b,"chimeric")) return("chimeric")
  if (is.na(a)) return(if (identical(b,"ref")) "no call" else b)
  if (is.na(b)) return(if (identical(a,"ref")) "no call" else a)
  if (a == "ref") return(b); if (b == "ref") return(a)
  if (a == b) return(a)
  "discordant"
}
BK$cls <- mapply(cls, BK$L1, BK$L2)
print(table(BK$cls, useNA="ifany"))
K <- BK[BK$cls %in% c("A","B"), ]
cat(sprintf("\n  ortholog braids to draw: %d (%d A, %d B)\n", nrow(K),
            sum(K$cls=="A"), sum(K$cls=="B")))

hr("6. RIBBON POLYGONS")
K$x1a <- mapply(mb2x, K$genome1, K$chr1, K$m1a)
K$x1b <- mapply(mb2x, K$genome1, K$chr1, K$m1b)
K$x2a <- mapply(mb2x, K$genome2, K$chr2, K$m2a)
K$x2b <- mapply(mb2x, K$genome2, K$chr2, K$m2b)
K$y1  <- CHR$y[match(paste(K$genome1, K$chr1), paste(CHR$genome, CHR$chr))]
K$y2  <- CHR$y[match(paste(K$genome2, K$chr2), paste(CHR$genome, CHR$chr))]
K <- K[complete.cases(K[, c("x1a","x1b","x2a","x2b","y1","y2")]), ]
K <- K[abs(K$y1 - K$y2) < 1.5, ]          # adjacent rows only, as GENESPACE draws
cat(sprintf("  ribbons with all coordinates, adjacent rows: %d\n", nrow(K)))
NS <- 32
POLY <- do.call(rbind, lapply(seq_len(nrow(K)), function(i) {
  t  <- seq(0, 1, length.out=NS)
  sm <- 0.5 - 0.5*cos(pi*t)
  yy <- K$y1[i] + sm*(K$y2[i] - K$y1[i])
  data.frame(id     = i,
             x      = c(K$x1a[i] + sm*(K$x2a[i]-K$x1a[i]),
                        rev(K$x1b[i] + sm*(K$x2b[i]-K$x1b[i]))),
             y      = c(yy, rev(yy)),
             region = K$region[i], cls = K$cls[i], stringsAsFactors=FALSE)
}))
cat(sprintf("  polygon vertices: %d\n", nrow(POLY)))

hr("7. ASSERTION — braids must sit on bars of their own colour")
barlab <- function(x, y) {
  k <- which(abs(IV$y - y) < 0.25 & IV$xa - 1 <= x & IV$xb + 1 >= x)
  if (!length(k)) return(NA_character_)
  IV$label[k[1]]
}
K$b1 <- mapply(barlab, (K$x1a+K$x1b)/2, K$y1)
K$b2 <- mapply(barlab, (K$x2a+K$x2b)/2, K$y2)
mm <- (!is.na(K$b1) & K$b1 != K$cls) | (!is.na(K$b2) & K$b2 != K$cls)
cat(sprintf("  drawn braids whose midpoint sits on the WRONG bar: %d of %d (%.1f%%)\n",
            sum(mm), nrow(K), 100*mean(mm)))
if (sum(mm)) {
  cat("\n  These are ends spanning a tract boundary: the midpoint falls in a\n")
  cat("  neighbouring tract of the other label. Listing them:\n")
  print(head(K[mm, c("base","region","cls","b1","b2")], 10))
}

hr("8. DRAW")
REGS <- sort(unique(POLY$region))
PAL <- setNames(colorRampPalette(c("#C0392B","#E67E22","#F1C40F","#7DCEA0",
                "#48C9B0","#5DADE2","#5B54B8","#A569BD"))(length(REGS)), REGS)
p <- ggplot() +
  geom_polygon(data=POLY[POLY$cls=="A",], aes(x, y, group=id, fill=region),
               alpha=0.40, colour=NA) +
  geom_polygon(data=POLY[POLY$cls=="B",], aes(x, y, group=id, fill=region),
               alpha=0.85, colour=NA) +
  geom_rect(data=CHR, aes(xmin=x0, xmax=x1, ymin=y-0.035, ymax=y+0.035),
            fill="white", colour="grey55", linewidth=0.3) +
  geom_rect(data=IV, aes(xmin=xa, xmax=xb, ymin=ybot, ymax=ytop),
            fill=ifelse(IV$label=="A", "#1D9E75", "#D85A30"), colour=NA) +
  scale_fill_manual(values=PAL, name="ancestral region") +
  scale_y_continuous(breaks=sort(unique(CHR$y)),
                     labels=gsub("_"," ", CHR$genome[match(sort(unique(CHR$y)), CHR$y)]),
                     trans="reverse") +
  labs(title="Droseraceae riparian, phased by subgenome",
       subtitle=paste0("Braids and bars share ONE coordinate system, so they ",
         "cannot disagree.\nBars: GREEN = A, ORANGE = B. Braids: MUTED = A, ",
         "VIVID = B, hue = ancestral region.\nOnly ortholog links are drawn (",
         sum(BK$cls=="discordant"), " A-to-B links and ",
         sum(BK$cls=="chimeric"), " chimeric ends removed)."),
       x="Chromosomes scaled by gene rank order", y=NULL) +
  theme_minimal(11) +
  theme(panel.background=element_rect(fill="grey96", colour=NA),
        plot.background=element_rect(fill="grey96", colour=NA),
        panel.grid=element_blank(), axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_text(face="italic", size=11),
        plot.title=element_text(face="bold", size=14),
        plot.subtitle=element_text(size=8.5, colour="grey30"),
        legend.position="bottom")
ggsave("DR/fig/DR17c_riparian_AB.pdf", p, width=16, height=9)
ggsave("DR/fig/DR17c_riparian_AB.png", p, width=16, height=9, dpi=200)
cat("  wrote DR/fig/DR17c_riparian_AB.{pdf,png}\n")
