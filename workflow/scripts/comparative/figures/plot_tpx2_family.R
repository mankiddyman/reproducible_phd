suppressMessages({library(ggtree); library(ape); library(ggplot2); library(data.table); library(Biostrings)})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/holocentricity/tpx2")
outdir <- "../figures"

tr <- read.tree("tpx2_family.nwk")
seqs <- readAAStringSet("tpx2_family.faa"); names(seqs) <- sub(" .*","",names(seqs))
meta <- data.table(label=tr$tip.label)
meta[, length := width(seqs)[match(label, names(seqs))]]
meta[, group := fcase(
  grepl("^paradoxa|^regia|^roseana", label), "holocentric",
  grepl("^aliciae|^binata|^capensis|^filiformis|^tokaiensis", label), "monocentric",
  grepl("^Atha", label), "Arabidopsis",
  grepl("^Beta", label), "Beta",
  grepl("^Dionaea", label), "Dionaea",
  grepl("^Nepenthes", label), "Nepenthes", default="other")]
meta[, cls := fifelse(length < 250, "short (<250 aa)", "long")]
cat("short seqs in the family:", meta[cls!="long", .N], "\n")
print(meta[cls!="long"][order(length), .(label, length, group)], nrows=40)

# the TPXL6-short clade = MRCA of the members we already verified
anchors <- intersect(c("paradoxa_g980.t1","regia_g4881.t1","roseana_h2tg000072l-g133.t1",
                       "Atha_AT5G37478","Nepenthes__BSYO01000005.1_000631.1",
                       "Dionaea__scaffold16_001058.1"), tr$tip.label)
cat("\nanchors found:", length(anchors), "of 6 ->", paste(anchors, collapse=", "), "\n")
n_short <- if(length(anchors) >= 2) getMRCA(tr, anchors) else NA
cat("short-lineage MRCA node:", n_short, "| tips in it:",
    if(!is.na(n_short)) length(extract.clade(tr, n_short)$tip.label) else NA, "\n")

# WVD2/WDL cortical-MT anchor, for orientation
n_wvd <- if("Atha_AT5G28646" %in% tr$tip.label) which(tr$tip.label=="Atha_AT5G28646") else NA

# --- which short seqs are IN the clade vs OUTSIDE it? ---
in_clade <- if(!is.na(n_short)) extract.clade(tr, n_short)$tip.label else character()
meta[, in_short_clade := label %in% in_clade]
cat("\n=== the", length(in_clade), "tips in the short clade ===\n")
print(meta[label %in% in_clade][order(length), .(label, length, group, cls)])
cat("\n=== SHORT seqs sitting OUTSIDE the clade (different short proteins) ===\n")
oo <- meta[cls!="long" & in_short_clade==FALSE][order(length), .(label, length, group)]
if(nrow(oo)) print(oo, nrows=40) else cat("  none - every short protein is in the clade\n")
cat("\n=== where is WVD2 (cortical MT)? ===\n")
if(!is.na(n_wvd)){
  sis <- extract.clade(tr, tr$edge[tr$edge[,2]==n_wvd, 1])$tip.label
  cat("  nearest neighbours:", paste(setdiff(sis, "Atha_AT5G28646"), collapse=", "), "\n")
  cat("  inside the short clade?", "Atha_AT5G28646" %in% in_clade, "\n")
} else cat("  AT5G28646 not in the tree\n")

cols <- c("holocentric"="#1b9e77","monocentric"="#d95f02","Arabidopsis"="#7570b3",
          "Beta"="#e7298a","Dionaea"="#66a61e","Nepenthes"="#e6ab02","other"="grey50")

p <- ggtree(tr, layout="circular", size=0.28) %<+% meta
if(!is.na(n_short)) p <- p + geom_hilight(node=n_short, fill="#fdbf6f", alpha=0.35)
p <- p +
  geom_tippoint(aes(color=group, shape=cls), size=1.5) +
  scale_shape_manual(values=c("long"=16, "short (<250 aa)"=17), name="protein length") +
  scale_color_manual(values=cols, name="lineage") +
  theme(legend.position="right", legend.text=element_text(size=9),
        plot.title=element_text(face="bold", size=16),
        plot.background=element_rect(fill="white", colour=NA)) +
  labs(title="The TPX2 domain family across Droseraceae and reference genomes",
       subtitle=sprintf("%d loci, 12 genomes. Highlighted = the short lineage containing TPXL6 (Arabidopsis AT5G37478).",
                        length(tr$tip.label)))
if(!is.na(n_short))
  p <- p + geom_cladelab(node=n_short, label="TPXL6-short", offset=0.15, fontsize=3.8, barsize=0.8)
if(!is.na(n_wvd)) p <- p + geom_cladelab(node=n_wvd, label="WVD2 (cortical MT)",
                                         offset=0.15, fontsize=3.2, barsize=0, textcolour="grey30")

ggsave(file.path(outdir,"tpx2_family_overview.pdf"), p, width=11, height=10)
ggsave(file.path(outdir,"tpx2_family_overview.jpg"), p, width=11, height=10, dpi=200)
cat("\nwrote tpx2_family_overview.{pdf,jpg}\n")
