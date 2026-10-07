#!/usr/bin/env Rscript
# ============================================================================
# DR17b — ZOOM ON ONE CHROMOSOME
#
# The tract table says capensis chr15_collapsed has chr6_dom = B over
# x 6200.2-6672.8 with no other tract in that window. The figure appears to
# show green inside it. This renders that one chromosome large enough to see
# exactly what is drawn where, with every tract boundary and every braid
# attachment point marked.
#
# If the render disagrees with the table, it is a drawing bug and it will be
# visible here.
#
# IN   ../genespace/riparian/Nepenthes_gracilis_geneOrder_rSourceData.rda
#      ../genespace/results/combBed.txt
#      DR/out/DR02_propagated_blocks.csv
# OUT  DR/fig/DR17b_zoom_capensis_chr15.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(GENESPACE); library(ggplot2)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
GS <- file.path(dirname(getwd()), "genespace")

GG <- "Drosera_capensis"; CC <- "chr15_collapsed"

e <- new.env()
load(file.path(GS,"riparian","Nepenthes_gracilis_geneOrder_rSourceData.rda"), envir=e)
P  <- e$srcd$ggplotObj
BR <- as.data.frame(P$layers[[1]]$data)
BR$base <- sub("_[0-9]+ chr[^ ]*$", "", BR$blkID)
CH <- as.data.frame(P$layers[[2]]$data)
CHR <- do.call(rbind, lapply(split(CH, paste(CH$genome,CH$chr)), function(z)
  data.frame(genome=z$genome[1], chr=z$chr[1], x0=min(z$x), x1=max(z$x),
             y=mean(z$y), stringsAsFactors=FALSE)))
BED <- read.table(file.path(GS,"results","combBed.txt"), header=TRUE, sep="\t",
                  quote="", comment.char="", stringsAsFactors=FALSE)
BED$mb <- (BED$start+BED$end)/2e6
G <- read.csv("DR/out/DR02_propagated_blocks.csv", stringsAsFactors=FALSE)
G$mb <- G$mid/1e6

k  <- which(CHR$genome==GG & CHR$chr==CC)
yv <- CHR$y[k]
bb <- BED[BED$genome==GG & BED$chr==CC, ]; bb <- bb[order(bb$ord), ]
mb2x <- function(mb) CHR$x0[k] +
  ((which.min(abs(bb$mb-mb))-1)/max(1,nrow(bb)-1))*(CHR$x1[k]-CHR$x0[k])

# tracts, exactly as DR17 builds them
z <- G[G$genome==GG & G$chr==CC, ]; z <- z[order(z$mb), ]
r <- rle(paste(z$region, z$label)); en <- cumsum(r$lengths); st <- en-r$lengths+1
TR <- do.call(rbind, lapply(seq_along(r$lengths), function(j){
  w <- z[st[j]:en[j], ]
  data.frame(region=w$region[1], label=w$label[1],
             lo=min(w$mb), hi=max(w$mb), n=nrow(w), stringsAsFactors=FALSE)}))
TR <- TR[TR$n >= 10, ]
TR$xa <- sapply(TR$lo, mb2x); TR$xb <- sapply(TR$hi, mb2x)
TR$col <- ifelse(TR$label=="A", "#1D9E75", "#D85A30")
cat("tracts drawn:\n"); print(TR[, c("region","label","lo","hi","xa","xb","n")])

# every braid attaching to this chromosome
ATT <- do.call(rbind, lapply(unique(BR$blkID), function(id) {
  q <- BR[BR$blkID==id, ]
  out <- NULL
  for (side in c("lo","hi")) {
    ee <- if (side=="lo") q[q$y < min(q$y)+0.02, ] else q[q$y > max(q$y)-0.02, ]
    yy <- if (side=="lo") min(q$y) else max(q$y)
    if (abs(yy-yv) > 0.25) next
    if (max(ee$x) < CHR$x0[k] || min(ee$x) > CHR$x1[k]) next
    out <- rbind(out, data.frame(blkID=id, region=sub("^.*_[0-9]+ ","",id),
                 xa=min(ee$x), xb=max(ee$x), stringsAsFactors=FALSE))
  }
  out
}))
cat(sprintf("\nbraid attachments: %d\n", nrow(ATT)))
ATT$y <- seq(-1.2, -3.6, length.out=nrow(ATT))
pal <- c(chr2_dom="#C0392B", chr3_dom="#E67E22", chr4_dom="#7DCEA0",
         chr5_dom="#48C9B0", chr6_dom="#5DADE2", chr7_dom="#5B54B8",
         chr9_dom="#A569BD", chr10_dom="#F1C40F")
ATT$col <- unname(pal[ATT$region]); ATT$col[is.na(ATT$col)] <- "grey50"

p <- ggplot() +
  geom_rect(aes(xmin=CHR$x0[k], xmax=CHR$x1[k], ymin=-0.35, ymax=0.35),
            fill="grey92", colour="grey60") +
  geom_rect(data=TR, aes(xmin=xa, xmax=xb, ymin=-0.3, ymax=0.3),
            fill=TR$col, colour="white", linewidth=0.3) +
  geom_text(data=TR, aes(x=(xa+xb)/2, y=0,
            label=paste0(sub("_dom","",region),"\n",label)),
            size=2.8, colour="white", fontface="bold") +
  geom_segment(data=ATT, aes(x=xa, xend=xb, y=y, yend=y),
               colour=ATT$col, linewidth=2.5) +
  geom_text(data=ATT, aes(x=(xa+xb)/2, y=y-0.28,
            label=sub("_dom","",region)), size=2.2, colour=ATT$col) +
  scale_x_continuous(breaks=pretty(c(CHR$x0[k], CHR$x1[k]), 12)) +
  labs(title=sprintf("%s %s — tracts (top) and braid attachments (below)", GG, CC),
       subtitle=paste0("Bar colour: GREEN = A, ORANGE = B, labelled with its ",
         "ancestral region.\nEach line below is one braid end, coloured by ",
         "its own ancestral region."),
       x="plot x coordinate", y=NULL) +
  theme_minimal(10) +
  theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        axis.text.y=element_blank())
ggsave("DR/fig/DR17b_zoom_capensis_chr15.pdf", p, width=16, height=7)
cat("\nwrote DR/fig/DR17b_zoom_capensis_chr15.pdf\n")
