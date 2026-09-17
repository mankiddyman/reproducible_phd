#!/usr/bin/env Rscript
# ============================================================================
# DR18 — RIPARIAN REBUILT FROM ALL BLOCKS, PHASED BY SUBGENOME
#
# WHY NOT REUSE GENESPACE'S SAVED PLOT
#   plot_riparian() calls phase_blks(), which selects a representative subset
#   for legibility: 1,278 braids out of 5,710 blocks between adjacent genomes.
#   Orthologous links that exist in the data were simply not drawn. This builds
#   every adjacent-pair block instead.
#
# WHY THE BARS ARE PER REGION
#   A Drosera chromosome is a FUSION of several ancestral regions, interleaved
#   along its length (57 of 68 chromosomes carry both A and B). Collapsing them
#   into one strip means the bar can read B from one region's genes while a
#   braid correctly reads A from another region's genes at the same x. Both
#   right, visually contradictory. So each chromosome gets one thin row PER
#   REGION, and a braid for region X sits directly under region X's row.
#
# THE FILTER
#   Every block belongs to exactly one ancestral region (the Nepenthes _dom
#   chromosome). A block is drawn only if BOTH ends resolve to the same
#   subgenome FOR THAT REGION. A-to-B links are homeologous, not orthologous,
#   and synteny cannot tell them apart in a polyploid.
#
# IN   ../genespace/results/{syntenicBlock_coordinates.csv,combBed.txt}
#      DR/out/DR02_propagated.csv
# OUT  DR/fig/DR18_riparian_full.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
GS <- file.path(dirname(getwd()), "genespace")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
ORD <- c("Nepenthes_gracilis","Dionaea_muscipula","Drosera_regia",
         "Drosera_binata","Drosera_paradoxa","Drosera_scorpioides",
         "Drosera_capensis")

hr("1. LAYOUT — chromosomes in gene rank order")
BED <- read.table(file.path(GS,"results","combBed.txt"), header=TRUE, sep="\t",
                  quote="", comment.char="", stringsAsFactors=FALSE)
BED <- BED[BED$genome %in% ORD, ]
BED <- BED[order(BED$genome, BED$chr, BED$ord), ]
CHR <- do.call(rbind, lapply(split(BED, paste(BED$genome, BED$chr)), function(z)
  data.frame(genome=z$genome[1], chr=z$chr[1], n=nrow(z),
             mbmin=min(z$start)/1e6, mbmax=max(z$end)/1e6,
             stringsAsFactors=FALSE)))
CHR <- CHR[CHR$n >= 100, ]
CHR <- CHR[order(CHR$genome, -CHR$n), ]
GAPP <- 0.012
CHR$x0 <- NA_real_; CHR$x1 <- NA_real_
for (g in ORD) {
  k <- which(CHR$genome==g); if (!length(k)) next
  tot <- sum(CHR$n[k]); gp <- GAPP*tot
  cum <- 0
  for (i in k) { CHR$x0[i] <- cum; CHR$x1[i] <- cum + CHR$n[i]; cum <- cum + CHR$n[i] + gp }
  w <- max(CHR$x1[k]); CHR$x0[k] <- CHR$x0[k]/w; CHR$x1[k] <- CHR$x1[k]/w
}
CHR$y <- match(CHR$genome, ORD)
cat(sprintf("  chromosomes drawn: %d\n", nrow(CHR)))
print(as.data.frame(table(CHR$genome)))

hr("2. GENE -> x, and the per-region label")
G <- read.csv("DR/out/DR02_propagated.csv", stringsAsFactors=FALSE)
G$mb <- G$mid/1e6
PAIRS <- read.csv("fractionation_by_chrpair.csv", stringsAsFactors=FALSE)
s1 <- PAIRS$retained_more == PAIRS$chrA
DIOSIDE <- c(setNames(ifelse(s1,"A","B"), PAIRS$chrA),
             setNames(ifelse(s1,"B","A"), PAIRS$chrB))
rank_of <- function(gg, cc, mb) {
  bb <- BED[BED$genome==gg & BED$chr==cc, ]
  if (!nrow(bb)) return(NA_real_)
  bb$rank <- seq_len(nrow(bb)); bb$mb <- (bb$start+bb$end)/2e6
  bb <- bb[order(bb$mb), ]
  bb$rank[pmax(1, findInterval(mb, bb$mb))] / nrow(bb)
}
tox <- function(gg, cc, mb) {
  k <- which(CHR$genome==gg & CHR$chr==cc); if (!length(k)) return(NA_real_)
  f <- rank_of(gg, cc, mb); if (is.na(f)) return(NA_real_)
  CHR$x0[k] + f * (CHR$x1[k] - CHR$x0[k])
}

hr("3. BLOCKS -> RIBBONS, filtered per region")
BC <- fread(file.path(GS,"results","syntenicBlock_coordinates.csv"))
BC <- as.data.frame(BC, stringsAsFactors=FALSE)

adj <- data.frame(a=ORD[-length(ORD)], b=ORD[-1], stringsAsFactors=FALSE)
keep <- rep(FALSE, nrow(BC))
for (i in seq_len(nrow(adj)))
  keep <- keep | (BC$genome1==adj$a[i] & BC$genome2==adj$b[i]) |
                 (BC$genome1==adj$b[i] & BC$genome2==adj$a[i])
BC <- BC[keep, ]
cat(sprintf("  adjacent-pair blocks with a region: %d\n", nrow(BC)))

# returns "region|label" for the genes inside a block interval
reg_lab <- function(gg, cc, m0, m1) {
  if (gg=="Nepenthes_gracilis") return(paste0(sub("_dom$","_dom",cc), "|ref"))
  if (gg=="Dionaea_muscipula")  return(paste0(NA, "|", unname(DIOSIDE[cc])))
  z <- G[G$genome==gg & G$chr==cc & G$mb>=m0 & G$mb<=m1, ]
  if (nrow(z) < 3) return(NA_character_)
  tr <- table(z$region); rg <- names(tr)[which.max(tr)]
  z <- z[z$region==rg, ]
  tl <- table(z$label)
  if (max(tl)/sum(tl) < 0.7) return(NA_character_)
  paste0(rg, "|", names(tl)[which.max(tl)])
}
m1a <- pmin(BC$startBp1,BC$endBp1)/1e6; m1b <- pmax(BC$startBp1,BC$endBp1)/1e6
m2a <- pmin(BC$startBp2,BC$endBp2)/1e6; m2b <- pmax(BC$startBp2,BC$endBp2)/1e6
r1 <- mapply(reg_lab, BC$genome1, BC$chr1, m1a, m1b)
r2 <- mapply(reg_lab, BC$genome2, BC$chr2, m2a, m2b)
BC$L1 <- sub("^.*\\|", "", r1); BC$L2 <- sub("^.*\\|", "", r2)
g1r <- sub("\\|.*$", "", r1);   g2r <- sub("\\|.*$", "", r2)
BC$region <- ifelse(!is.na(g1r) & g1r != "NA", g1r, g2r)
BC$L1[is.na(r1)] <- NA_character_; BC$L2[is.na(r2)] <- NA_character_
cat(sprintf("  blocks with a region resolved: %d of %d\n",
            sum(!is.na(BC$region) & BC$region != "NA"), nrow(BC)))
# both ends must agree on the ancestral region, or they are not comparable
mismatch <- !is.na(g1r) & !is.na(g2r) & g1r != "NA" & g2r != "NA" & g1r != g2r
cat(sprintf("  ends disagreeing on region: %d (dropped)\n", sum(mismatch)))
BC$L1[mismatch] <- NA_character_
cls <- ifelse(BC$L1=="ref" & !is.na(BC$L2), BC$L2,
       ifelse(BC$L2=="ref" & !is.na(BC$L1), BC$L1,
       ifelse(is.na(BC$L1)|is.na(BC$L2), NA_character_,
       ifelse(BC$L1==BC$L2, BC$L1, "discordant"))))
BC$cls <- cls
BC$m1a <- m1a; BC$m1b <- m1b; BC$m2a <- m2a; BC$m2b <- m2b
cat("\n  blocks by class:\n"); print(table(BC$cls, useNA="ifany"))
KEEP <- BC[!is.na(BC$cls) & BC$cls %in% c("A","B"), ]
cat(sprintf("\n  ORTHOLOG blocks drawn: %d (%d A, %d B)\n", nrow(KEEP),
            sum(KEEP$cls=="A"), sum(KEEP$cls=="B")))
cat(sprintf("  homeolog links dropped: %d\n", sum(BC$cls=="discordant", na.rm=TRUE)))

hr("4. BUILD THE POLYGONS")
KEEP$xa0 <- mapply(tox, KEEP$genome1, KEEP$chr1, KEEP$m1a)
KEEP$xa1 <- mapply(tox, KEEP$genome1, KEEP$chr1, KEEP$m1b)
KEEP$xb0 <- mapply(tox, KEEP$genome2, KEEP$chr2, KEEP$m2a)
KEEP$xb1 <- mapply(tox, KEEP$genome2, KEEP$chr2, KEEP$m2b)
KEEP$y1 <- match(KEEP$genome1, ORD); KEEP$y2 <- match(KEEP$genome2, ORD)
KEEP <- KEEP[complete.cases(KEEP[,c("xa0","xa1","xb0","xb1","y1","y2")]), ]
cat(sprintf("  ribbons with all coordinates: %d\n", nrow(KEEP)))
NS <- 30
POLY <- do.call(rbind, lapply(seq_len(nrow(KEEP)), function(i) {
  t <- seq(0,1,length.out=NS); sm <- 0.5-0.5*cos(pi*t)
  y1 <- KEEP$y1[i]; y2 <- KEEP$y2[i]
  data.frame(id=i,
    x=c(KEEP$xa0[i]+sm*(KEEP$xb0[i]-KEEP$xa0[i]),
        rev(KEEP$xa1[i]+sm*(KEEP$xb1[i]-KEEP$xa1[i]))),
    y=c(y1+sm*(y2-y1), rev(y1+sm*(y2-y1))),
    region=KEEP$region[i], cls=KEEP$cls[i], stringsAsFactors=FALSE)
}))
cat(sprintf("  polygon vertices: %d\n", nrow(POLY)))

hr("5. PER-REGION BARS")
BARS <- do.call(rbind, lapply(split(G, paste(G$genome,G$chr,G$region)), function(z){
  if (nrow(z) < 10) return(NULL)
  z <- z[order(z$mb), ]
  do.call(rbind, lapply(split(z, c(0,cumsum(diff(z$mb) > 0.5))), function(w){
    if (nrow(w) < 8) return(NULL)
    data.frame(genome=w$genome[1], chr=w$chr[1], region=w$region[1],
               label=w$label[1], m0=min(w$mb), m1=max(w$mb), stringsAsFactors=FALSE)
  }))
}))
BARS$xa <- mapply(tox, BARS$genome, BARS$chr, BARS$m0)
BARS$xb <- mapply(tox, BARS$genome, BARS$chr, BARS$m1)
BARS <- BARS[!is.na(BARS$xa) & !is.na(BARS$xb) & BARS$label %in% c("A","B"), ]
REGS <- sort(unique(G$region))
BARS$slot <- match(BARS$region, REGS)
BARS$y <- match(BARS$genome, ORD) - 0.06 - 0.028*BARS$slot
DB <- CHR[CHR$genome=="Dionaea_muscipula", ]
DB$label <- unname(DIOSIDE[DB$chr]); DB <- DB[!is.na(DB$label), ]
cat(sprintf("  region bars: %d | Dionaea chromosome bars: %d\n", nrow(BARS), nrow(DB)))

hr("6. DRAW")
PAL <- setNames(colorRampPalette(c("#C0392B","#E67E22","#F1C40F","#7DCEA0",
                                   "#48C9B0","#5DADE2","#5B54B8","#A569BD"))(length(REGS)), REGS)
p <- ggplot() +
  geom_polygon(data=POLY[POLY$cls=="A",],
               aes(x,y,group=id,fill=region), alpha=0.35, colour=NA) +
  geom_polygon(data=POLY[POLY$cls=="B",],
               aes(x,y,group=id,fill=region), alpha=0.80, colour=NA) +
  geom_rect(data=CHR, aes(xmin=x0,xmax=x1,ymin=y-0.035,ymax=y+0.035),
            fill="grey25", colour=NA) +
  geom_rect(data=BARS, aes(xmin=xa,xmax=xb,ymin=y-0.011,ymax=y+0.011,
                           colour=label), fill=NA, linewidth=0.9) +
  geom_rect(data=DB, aes(xmin=x0,xmax=x1,ymin=y-0.075,ymax=y-0.045,
                         colour=label), fill=NA, linewidth=1.1) +
  scale_fill_manual(values=PAL, name="ancestral region") +
  scale_colour_manual(values=c(A="#1D9E75", B="#D85A30"), name="subgenome") +
  scale_y_continuous(breaks=seq_along(ORD), labels=gsub("_"," ",ORD),
                     trans="reverse") +
  labs(title="Droseraceae riparian, all blocks, phased by subgenome",
       subtitle=paste0(nrow(KEEP), " ortholog blocks drawn from ", nrow(BC),
         " adjacent-pair blocks (GENESPACE's own riparian draws 1,278).\n",
         sum(BC$cls=="discordant", na.rm=TRUE),
         " A-to-B links removed as homeologous. MUTED = A, VIVID = B.\n",
         "Bars under each chromosome: one thin row PER ancestral region, ",
         "green = A, orange = B."),
       x="Chromosomes scaled by gene rank order", y=NULL) +
  theme_minimal(11) +
  theme(panel.grid=element_blank(), axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_text(face="italic", size=11),
        plot.title=element_text(face="bold", size=14),
        plot.subtitle=element_text(size=8.5, colour="grey30"),
        legend.position="bottom")
ggsave("DR/fig/DR18_riparian_full.pdf", p, width=17, height=11)
ggsave("DR/fig/DR18_riparian_full.png", p, width=17, height=11, dpi=200)
cat("  wrote DR/fig/DR18_riparian_full.{pdf,png}\n")
