suppressMessages({library(ggplot2); library(data.table)})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/holocentricity")
outdir <- "figures"; QLEN <- 140
GRN <- "#1b9e77"; ORG <- "#d95f02"

union_len <- function(qs, qe){
  if(!length(qs)) return(0L)
  o <- order(qs); qs <- qs[o]; qe <- qe[o]
  tot <- 0L; cs <- qs[1]; ce <- qe[1]
  if(length(qs) > 1) for(i in 2:length(qs)){
    if(qs[i] <= ce+1) ce <- max(ce, qe[i]) else { tot <- tot + (ce-cs+1); cs <- qs[i]; ce <- qe[i] }
  }
  as.integer(tot + (ce-cs+1))
}

# ---- genome: cluster HSPs into loci (same contig, <50 kb apart), take the best locus ----
g <- fread("tpx2/tpxl6_genome_hsps.tsv",
           col.names=c("sp","ctg","qstart","qend","sstart","send","frame","evalue","pident"))
g[, gmin := pmin(sstart, send)]
setorder(g, sp, ctg, gmin)
g[, newloc := as.integer(ctg != shift(ctg, fill="") | (gmin - shift(gmin, fill=-1e9)) > 50000), by=sp]
g[, locus := cumsum(newloc), by=sp]
loc <- g[, .(aa=union_len(qstart,qend), nframe=uniqueN(frame), nhsp=.N,
             best_e=min(evalue), pid=round(max(pident)), ctg=ctg[1]), by=.(sp,locus)]
setorder(loc, sp, -aa)
best <- loc[, .SD[1], by=sp]

# ---- annotation: union over short annotated hits ----
pa <- tryCatch({
  p <- fread("tpx2/tpxl6_prot_hits.tsv",
             col.names=c("sp","sseqid","qstart","qend","slen","evalue","pident"))
  p[, .(aa=union_len(qstart,qend)), by=sp]
}, error=function(e) data.table(sp=character(), aa=integer()))

info <- data.table(
  sp   = c("regia","paradoxa","roseana","binata","capensis","aliciae","tokaiensis","filiformis"),
  nice = c("D. regia","D. paradoxa","D. roseana","D. binata",
           "D. capensis","D. aliciae","D. tokaiensis","D. filiformis"),
  centromere = c("holocentric","holocentric","holocentric","monocentric",
                 "monocentric","monocentric","monocentric","monocentric"))
d <- merge(info, best[, .(sp, g_aa=aa, nframe, pid, ctg)], by="sp", all.x=TRUE, sort=FALSE)
d <- merge(d, pa[, .(sp, p_aa=aa)], by="sp", all.x=TRUE, sort=FALSE)
d[is.na(g_aa), g_aa := 0L]; d[is.na(p_aa), p_aa := 0L]
d[, y := .N:1]
cat("\n=== per species ===\n"); print(d[, .(nice, p_aa, g_aa, nframe, pid, ctg)])

XT <- 168
d[, annot := fifelse(p_aa > 0, "YES", "no")]
d[, note := fcase(
  sp %in% c("regia","paradoxa","roseana"),
    "intact protein-coding gene",
  sp == "binata",
    sprintf("as much DNA as the holocentrics \u2014 but split across %d\nreading frames with stop codons throughout.\nNo protein can be made from it.", nframe),
  sp == "capensis",
    "only a short domain shared with a different,\nmuch longer gene. The gene itself is gone.",
  default = "")]

p <- ggplot(d, aes(y=y)) +
  annotate("rect", xmin=-4, xmax=250, ymin=4.5, ymax=5.5, fill="grey94") +
  geom_segment(aes(x=0, xend=g_aa, yend=y, colour=centromere), linewidth=1.4, alpha=0.35) +
  geom_point(aes(x=g_aa, colour=centromere, fill=centromere,
                 shape=annot), size=5, stroke=1.4) +
  scale_shape_manual(values=c(YES=21, no=4), guide="none") +
  geom_text(aes(x=g_aa+4, label=paste0(g_aa, " aa")), hjust=0, size=3.4, colour="grey35") +
  geom_text(aes(x=XT, label=nice), hjust=0, size=4, fontface="italic",
            colour=ifelse(d$centromere=="holocentric", GRN, ORG)) +
  geom_text(aes(x=XT+36, label=annot, fontface=fifelse(annot=="YES","plain","bold")),
            hjust=0.5, size=3.8, colour=ifelse(d$annot=="YES","grey35",ORG)) +
  geom_text(aes(x=XT+52, label=note), hjust=0, size=3.05, colour="grey40", lineheight=0.95) +
  annotate("text", x=XT+36, y=8.95, label="annotated\nas a gene?", size=3.1,
           fontface="bold", colour="grey25", lineheight=0.9) +
  scale_colour_manual(values=c(holocentric=GRN, monocentric=ORG), name=NULL) +
  scale_fill_manual(values=c(holocentric=GRN, monocentric=ORG), guide="none") +
  scale_x_continuous("amino acids of the D. regia TPXL6 protein found by searching the raw DNA",
                     breaks=seq(0, 140, 35), limits=c(-6, 400), expand=c(0,0)) +
  scale_y_continuous(NULL, breaks=NULL, limits=c(0.4, 9.4), expand=c(0,0)) +
  labs(title="In D. binata the DNA is still there \u2014 but it is no longer a gene",
       subtitle="Every genome searched directly, ignoring the gene annotations. D. binata carries as much of the sequence as the\nholocentrics do \u2014 yet it is the only genome with that much sequence and no gene.",
       caption=paste0("Query: the intact D. regia TPXL6 protein (140 aa), tblastn against each genome; best locus per species (HSPs merged, overlaps counted once).\n",
                      "Even intact genes fall short of 140 aa because tblastn breaks at introns \u2014 so the holocentric values are the ceiling here, not 140.\n",
                      "Annotated = a short (<250 aa) protein in the species' annotation matches the query; the long TPX2-family paralogs are a different gene.\n",
                      "D. binata relic: 83% identical, present on both haplotypes \u2014 not an assembly artifact.")) +
  theme_minimal(base_size=13) +
  theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        axis.title.x=element_text(size=11, colour="grey30", margin=margin(t=8)),
        legend.position="none",
        plot.title=element_text(face="bold", size=17),
        plot.subtitle=element_text(size=10.5, colour="grey35", lineheight=1.1, margin=margin(b=20)),
        plot.caption=element_text(size=7.2, colour="grey45", hjust=0, lineheight=1.2),
        plot.background=element_rect(fill="white", colour=NA),
        plot.margin=margin(12,14,10,12))

ggsave(file.path(outdir,"tpxl6_dumbbell.pdf"), p, width=13, height=6.2)
ggsave(file.path(outdir,"tpxl6_dumbbell.jpg"), p, width=13, height=6.2, dpi=200)
cat("\nwrote tpxl6_dumbbell.{pdf,jpg}\n")
