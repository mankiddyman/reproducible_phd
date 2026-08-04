suppressMessages({library(ggplot2); library(data.table)})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/holocentricity")
outdir <- "figures"
GRN <- "#1b9e77"; ORG <- "#d95f02"

d <- fread("tpx2/tpxl6_genome_recovery.tsv")
print(d)

info <- data.table(
  sp        = c("regia","paradoxa","roseana","binata","capensis","aliciae","tokaiensis","filiformis"),
  nice      = c("D. regia","D. paradoxa","D. roseana","D. binata",
                "D. capensis","D. aliciae","D. tokaiensis","D. filiformis"),
  centromere= c("holocentric","holocentric","holocentric","monocentric",
                "monocentric","monocentric","monocentric","monocentric"),
  annotated = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  verdict   = c("intact gene", "intact gene", "intact gene",
                "PSEUDOGENE", "gone", "gone", "gone", "gone"))
d <- merge(info, d, by="sp", sort=FALSE)
d[, y := .N:1]
d[, is_pseudo := verdict=="PSEUDOGENE"]

# right-hand annotation, plain language
d[, note := fcase(
  verdict=="intact gene", sprintf("intact gene \u00b7 one reading frame \u00b7 annotated"),
  verdict=="PSEUDOGENE",  sprintf("just as much sequence as the holocentrics \u2014\nbut split across %d reading frames, stops throughout,\nno protein. Not annotated because there is no gene.", nframe),
  default =               "nothing but the shared domain of a\ndifferent, longer gene")]
d[verdict=="gone" & sp!="capensis", note := ""]     # write it once for the block

XT <- 86     # where the notes start

p <- ggplot(d, aes(x=qcov, y=y)) +
  annotate("rect", xmin=0, xmax=45, ymin=0.35, ymax=8.65, fill="grey93") +
  annotate("text", x=22.5, y=8.9, label="background: a different gene", size=3.1, color="grey45") +
  geom_col(aes(fill=centromere, linetype=is_pseudo), width=0.62,
           color="grey20", linewidth=0.5, orientation="y") +
  geom_text(aes(x=qcov-1.5, label=paste0(qcov, "%")),
            hjust=1, size=4.2, fontface="bold", color="white") +
  geom_text(aes(x=XT, label=verdict, fontface=fifelse(is_pseudo,"bold","plain")),
            hjust=0, size=4, color="grey10") +
  geom_text(aes(x=XT+26, y=y, label=note), hjust=0, size=3.05,
            color="grey40", lineheight=0.95) +
  scale_fill_manual(values=c(holocentric=GRN, monocentric=ORG), name=NULL) +
  scale_linetype_manual(values=c(`FALSE`="solid", `TRUE`="21"), guide="none") +
  scale_x_continuous("how much of the D. regia TPXL6 protein is recoverable from each genome",
                     breaks=c(0,25,50,75,100), labels=function(x) paste0(x,"%"),
                     limits=c(0,205), expand=c(0,0)) +
  scale_y_continuous(NULL, breaks=d$y, labels=d$nice, limits=c(0.35, 9.2), expand=c(0,0)) +
  labs(title="In D. binata the sequence is still there \u2014 but it is no longer a gene",
       subtitle="Searching each genome directly (tblastn), not the gene annotation. One query throughout: the intact D. regia protein.",
       caption=paste0("Query: D. regia g4881, 140 aa. Bars = % of it recovered from the best-matching contig. ",
                      "Holocentric bars fall short of 100% because tblastn breaks at introns.\n",
                      "D. binata relic: chr16, 83% identical, mirrored on both haplotypes \u2014 not an assembly artifact. ",
                      "The four 'gone' species hit only annotated long paralogs.\n",
                      "Nepenthes, Dionaea and Arabidopsis all carry an intact short copy. Two independent losses \u2014 n=2, no statistics.")) +
  theme_minimal(base_size=13) +
  theme(axis.text.y=element_text(face="italic", size=12),
        axis.title.x=element_text(size=11, color="grey30", margin=margin(t=8)),
        panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        legend.position=c(0.62, 1.06), legend.direction="horizontal",
        plot.title=element_text(face="bold", size=17),
        plot.subtitle=element_text(size=10.5, color="grey35", margin=margin(b=16)),
        plot.caption=element_text(size=7.4, color="grey45", hjust=0, lineheight=1.2),
        plot.background=element_rect(fill="white", color=NA),
        plot.margin=margin(12,14,10,12))

ggsave(file.path(outdir,"tpxl6_recovery.pdf"), p, width=13, height=6)
ggsave(file.path(outdir,"tpxl6_recovery.jpg"), p, width=13, height=6, dpi=200)
cat("wrote tpxl6_recovery.{pdf,jpg}\n")
