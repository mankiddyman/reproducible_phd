suppressMessages({library(ggplot2); library(data.table)})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/holocentricity")
outdir <- "figures"
GRN <- "#1b9e77"; ORG <- "#d95f02"
WHITE <- theme(plot.background=element_rect(fill="white", color=NA))

sp_lv <- rev(c("Nepenthes gracilis","Dionaea muscipula","D. regia","D. binata",
               "D. paradoxa","D. roseana","D. capensis","D. aliciae",
               "D. tokaiensis","D. filiformis"))
cent <- c("Nepenthes gracilis"="monocentric","Dionaea muscipula"="monocentric",
          "D. regia"="holocentric","D. binata"="monocentric","D. paradoxa"="holocentric",
          "D. roseana"="holocentric","D. capensis"="monocentric","D. aliciae"="monocentric",
          "D. tokaiensis"="monocentric","D. filiformis"="monocentric")
ypos <- setNames(seq_along(sp_lv), sp_lv)          # numeric y, no discrete scale

# ---------------- A: KNL2 (SANTA domain, tblastn vs genomes) ----------------
a <- data.table(sp = names(cent),
  hits  = c(117, 18, 14, 168, 50, 117, 23, 15, 12, 15),
  santa = c( 22,  0,  0,   0,  0,   0,  0,  0,  0,  0))
a[, `:=`(y = ypos[sp], centromere = cent[sp])]

pA <- ggplot(a, aes(x=santa, y=y)) +
  geom_col(aes(fill=centromere), width=0.65, orientation="y") +
  geom_text(aes(x=0.4, label=sprintf("%d tblastn hits searched \u2192 %d with SANTA", hits, santa)),
            hjust=0, size=3.2, color="grey25") +
  scale_fill_manual(values=c(holocentric=GRN, monocentric=ORG), guide="none") +
  scale_x_continuous("SANTA domains (PF09133) recovered", limits=c(0, 26), expand=c(0,0)) +
  scale_y_continuous(NULL, breaks=ypos, labels=names(ypos), limits=c(0.4, 10.6), expand=c(0,0)) +
  labs(subtitle="A.  \u03b1/\u03b2KNL2 \u2014 the search ran in every genome; only Nepenthes carries the domain") +
  theme_minimal(base_size=12) +
  theme(axis.text.y=element_text(face="italic", size=10.5),
        panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        plot.subtitle=element_text(face="bold", size=11.5)) + WHITE

# ---------------- B: TPXL6-short, proteome vs genome ----------------
b <- rbind(
  data.table(sp=names(cent), evidence="proteome (blastp, annotation-dependent)",
             nlp   = c(31, 31, 47, 0, 40, 43, 0, 0.39, 0, 2.2),
             nohit = c( F,  F,  F, T,  F,  F, T,    F, F,   F)),
  data.table(sp=names(cent), evidence="genome (tblastn, annotation-independent)",
             nlp   = c(NA, NA, 30, 23, 24, 22, 8, 7, 7, 6),
             nohit = FALSE))
b[, `:=`(y = ypos[sp], centromere = cent[sp])]
b <- b[!is.na(nlp)]

pB <- ggplot(b, aes(x=nlp, y=y)) +
  annotate("rect", xmin=-0.6, xmax=10, ymin=0.4, ymax=10.6, fill="grey92") +
  annotate("text", x=4.7, y=10.2, size=2.9, color="grey45",
           label="noise / cross-hits to long paralogs") +
  geom_line(aes(group=y), color="grey70", linewidth=0.4, orientation="y") +
  geom_point(aes(color=centromere, shape=evidence), size=3.6) +
  geom_text(data=b[nohit==TRUE], aes(x=0, label="no hit"),
            vjust=-1.2, size=2.7, color="grey40") +
  annotate("segment", x=23, xend=23, y=6.5, yend=6.85, colour="grey30",
           arrow=arrow(length=unit(0.18,"cm"))) +
  annotate("label", x=23, y=6.25, size=3.1, colour="grey15", fill="white",
           label="D. binata: absent from the proteome,\npresent in the genome = the pseudogene") +
  scale_color_manual(values=c(holocentric=GRN, monocentric=ORG), name="centromere") +
  scale_shape_manual(values=c(16, 17), name="evidence") +
  scale_x_continuous(expression(-log[10]~"(best e-value)"), limits=c(-0.6, 50), expand=c(0,0)) +
  scale_y_continuous(NULL, breaks=ypos, labels=names(ypos), limits=c(0.4, 10.6), expand=c(0,0)) +
  labs(subtitle="B.  TPXL6-short \u2014 the two evidence types agree everywhere except D. binata") +
  theme_minimal(base_size=12) +
  theme(axis.text.y=element_text(face="italic", size=10.5),
        panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        legend.position="bottom", legend.box="horizontal",
        plot.subtitle=element_text(face="bold", size=11.5)) + WHITE

CAP <- paste0("A: KNL2 \u03b1+\u03b2 queries, tblastn vs genomes (evalue 10), top-30 hit windows translated in 6 frames, hmmsearch PF09133. Nepenthes = positive control.\n",
              "B: AT5G37478 (178 aa) vs proteomes, all ranks, hits <250 aa; and holo proteins + AT5G37478 vs genomes. Outgroup genomes not tblastn-tested for TPXL6.\n",
              "Values <1 clamped to 0. Beta vulgaris: no short homolog (best 0.34, proteome only \u2014 not genome-verified).")

if (requireNamespace("patchwork", quietly=TRUE)) {
  suppressMessages(library(patchwork))
  p <- pA / pB + plot_layout(heights=c(1, 1.15)) +
    plot_annotation(title="Absence is only evidence when the positive control fires", caption=CAP,
      theme=theme(plot.title=element_text(face="bold", size=17),
                  plot.caption=element_text(size=7.3, colour="grey45", hjust=0, lineheight=1.15),
                  plot.background=element_rect(fill="white", color=NA)))
  ggsave(file.path(outdir,"detection_asymmetry.pdf"), p, width=12, height=8.5)
  ggsave(file.path(outdir,"detection_asymmetry.jpg"), p, width=12, height=8.5, dpi=200)
  cat("wrote detection_asymmetry.{pdf,jpg}\n")
} else {
  ggsave(file.path(outdir,"detection_asymmetry_A.jpg"), pA, width=12, height=4, dpi=200)
  ggsave(file.path(outdir,"detection_asymmetry_B.jpg"), pB, width=12, height=4.6, dpi=200)
  cat("patchwork missing -> wrote _A and _B separately\n")
}
