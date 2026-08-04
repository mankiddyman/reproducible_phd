suppressMessages({library(ggtree); library(ape); library(ggplot2); library(data.table); library(Biostrings)})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/holocentricity/tpx2")
outdir <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/holocentricity/figures"

tr <- read.tree("tpx2_short.treefile")

# root on the LONG-paralog clade (100 support): long genes + the two degraded holo long members
root_tips <- c("Atha_AT5G15510","regia_g4900.t1","capensis_g4300.t1",
               "paradoxa_g28615.t1","roseana_h1tg000004l-g260.t1")
root_tips <- intersect(root_tips, tr$tip.label)
cat("rooting on:", paste(root_tips, collapse=", "), "\n")
tr <- root(tr, outgroup=root_tips, resolve.root=TRUE)

seqs <- readAAStringSet("tpx2_short.faa"); names(seqs) <- sub(" .*","",names(seqs))
meta <- data.table(label=tr$tip.label, length=width(seqs)[match(tr$tip.label, names(seqs))])
meta[, group := fcase(
  grepl("^paradoxa|^regia|^roseana", label), "holocentric",
  grepl("^binata|^capensis", label), "monocentric",
  grepl("^Atha", label), "Arabidopsis",
  grepl("^Nepenthes", label), "Nepenthes",
  grepl("^Dionaea", label), "Dionaea", default="other")]
meta[, pseudo := grepl("PSEUDOGENE", label)]   # = the non-expressed binata copy
meta[, lineage := fifelse(label %in% root_tips, "long paralog", "short lineage")]
meta[, nice := fcase(
  grepl("PSEUDOGENE", label), "D. binata (not expressed)",
  grepl("^paradoxa", label), sub("paradoxa_","D. paradoxa ", sub("\\.t1$","",label)),
  grepl("^regia",    label), sub("regia_","D. regia ", sub("\\.t1$","",label)),
  grepl("^roseana",  label), sub("roseana_[^-]*-","D. roseana ", sub("\\.t1$","",label)),
  grepl("^capensis", label), sub("capensis_","D. capensis ", sub("\\.t1$","",label)),
  label=="Atha_AT5G37478", "Arabidopsis AT5G37478",
  label=="Atha_AT5G15510", "Arabidopsis AT5G15510",
  grepl("^Nepenthes", label), "Nepenthes gracilis",
  grepl("^Dionaea",   label), "Dionaea muscipula", default=label)]
meta[, tiplab := paste0(nice, "  (", length, "aa)")]

cols <- c("holocentric"="#1b9e77","monocentric"="#d95f02","Arabidopsis"="#7570b3",
          "Nepenthes"="#e6ab02","Dionaea"="#66a61e","other"="grey50")
maxd <- max(node.depth.edgelength(tr))

p <- ggtree(tr, size=0.7) %<+% meta +
  geom_tippoint(aes(color=group, shape=pseudo), size=3.4) +
  scale_shape_manual(values=c(`FALSE`=16, `TRUE`=8), guide="none") +
  geom_tiplab(aes(color=group, label=tiplab,
                  fontface=fifelse(pseudo,"bold.italic","italic")),
              size=3.8, offset=maxd*0.02) +
  geom_text2(aes(subset=!isTip & !is.na(suppressWarnings(as.numeric(label))),
                 label=label), hjust=1.25, vjust=-0.5, size=2.7, color="grey45") +
  scale_color_manual(values=cols, name=NULL) +
  xlim(0, maxd*2.0) +
  labs(title="Short TPX2-domain lineage, rooted on long paralogs",
       subtitle="Node labels = ultrafast bootstrap (%). 123 alignment columns \u2014 deep nodes poorly supported. \u2733 = intact ORF, not detected in leaf RNA.") +
  theme_tree() +
  theme(legend.position=c(0.12,0.78), legend.text=element_text(size=11),
        plot.title=element_text(face="bold", size=16),
        plot.subtitle=element_text(size=10, color="grey35"))

ggsave(file.path(outdir,"tpx2_short_tree.pdf"), p, width=12, height=7)
ggsave(file.path(outdir,"tpx2_short_tree.jpg"), p, width=12, height=7, dpi=200)
cat("wrote tpx2_short_tree.{pdf,jpg}\n")
print(meta[, .(nice, group, length, lineage)])
