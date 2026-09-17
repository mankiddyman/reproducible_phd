#!/usr/bin/env Rscript
# ============================================================================
# DR05h — THE 13-TIP TREE, real branch lengths, fully annotated
#
# One picture of the whole model:
#   - A tips and B tips form TWO CLADES  -> the allopolyploidy itself
#   - within each, the species relationships
#   - regia's position in each subgenome, and where they disagree
#
# Branch lengths are substitutions per site from the concatenated ML fit
# (GTR+F+I+G4, 2,064,141 sites, 1232 loci). Unlike ASTRAL's coalescent units
# these are real and comparable across the tree -- which is why this figure
# uses concatenation despite ASTRAL being the better topology estimate at
# regia's node. Weak nodes are flagged by sCF so the reader can see which
# splits the branch lengths are hanging off.
#
# ANNOTATION
#   UFBoot/SH-aLRT from the treefile; sCF from the --scf run. sCF is the
#   honest number: bootstrap saturates at 100 on 2 Mb of concatenated data,
#   sCF says what fraction of SITES actually support each branch. Nodes below
#   40% sCF are drawn in red.
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(ape); library(ggtree); library(ggplot2); library(dplyr); library(tidyr)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))

t <- read.tree("DR/tree/concat/all13_full.treefile")
cat("tips:", length(t$tip.label), "\n")

# --- sCF: match by descendant tip set, which survives rerooting ---
cf <- read.table("DR/tree/concat/all13_full_scf.cf.stat", header=TRUE, comment.char="#")
tipset <- function(tr, node) sort(tr$tip.label[Descendants(tr, node)])
Descendants <- function(tr, node) {
  if (node <= length(tr$tip.label)) return(node)
  ch <- tr$edge[tr$edge[,1]==node, 2]
  unlist(lapply(ch, function(x) Descendants(tr, x)))
}
pre <- lapply(cf$ID, function(n) if (n > length(t$tip.label)) tipset(t, n) else NA)
names(pre) <- cf$ID

tr <- root(t, outgroup="Nepenthes", resolve.root=TRUE)
tr$edge.length[is.na(tr$edge.length)] <- 0
nodes <- (length(tr$tip.label)+1):(length(tr$tip.label)+tr$Nnode)
post <- lapply(nodes, function(n) tipset(tr, n)); names(post) <- nodes

ann <- bind_rows(lapply(seq_along(post), function(i) {
  s <- post[[i]]; alt <- setdiff(tr$tip.label, s)
  hit <- which(vapply(pre, function(p)
    !all(is.na(p)) && (identical(sort(p), s) || identical(sort(p), sort(alt))),
    logical(1)))
  if (!length(hit)) return(NULL)
  tibble(node=as.integer(names(post)[i]), sCF=cf$sCF[cf$ID==as.integer(names(pre)[hit[1]])])
}))
cat(sprintf("sCF matched to %d of %d internal nodes\n", nrow(ann), length(nodes)))
print(as.data.frame(ann %>% arrange(sCF)), row.names=FALSE)

# --- support values from the treefile node labels ---
sup <- tibble(node = nodes,
              lab  = tr$node.label[nodes - length(tr$tip.label)]) %>%
       mutate(lab = ifelse(is.na(lab) | lab=="" | lab=="Root", NA, lab))

D <- tibble(label = tr$tip.label) %>%
  mutate(sub = ifelse(grepl("_A$", label), "A subgenome",
               ifelse(grepl("_B$", label), "B subgenome", "outgroup")),
         sp  = sub("_[AB]$", "", label),
         pretty = ifelse(sp=="Nepenthes", "Nepenthes gracilis",
                  ifelse(sp=="Dionaea", "Dionaea muscipula",
                         paste0("D. ", sp))),
         pretty = ifelse(sub=="outgroup", pretty,
                         paste0(pretty, "  ", substr(sub,1,1))))

NODEDAT <- full_join(sup, ann, by="node") %>%
  mutate(weak = !is.na(sCF) & sCF < 40,
         txt  = ifelse(is.na(sCF), lab, sprintf("%s | sCF %.0f", lab, sCF)))

p <- ggtree(tr, size=0.75) %<+% D +
  geom_tiplab(aes(label=pretty, colour=sub), fontface=3, size=3.6, offset=0.008) +
  scale_colour_manual(name=NULL,
    values=c("A subgenome"="#1D9E75", "B subgenome"="#D85A30",
             "outgroup"="grey40")) +
  ggnewscale::new_scale_colour()
p <- p %<+% NODEDAT +
  geom_nodepoint(aes(colour=weak), size=2) +
  geom_text2(aes(subset=!isTip, label=txt, colour=weak),
             hjust=1.08, vjust=-0.6, size=2.4) +
  scale_colour_manual(values=c("FALSE"="grey40","TRUE"="#B03020"), guide="none") +
  geom_treescale(x=0, y=-0.6, width=0.05, fontsize=3, linesize=0.5,
                 offset=0.15) +
  labs(title = "Droseraceae subgenome phylogeny",
       subtitle = paste0(
         "Concatenated ML, GTR+F+I+G4, 1,232 loci / 2,064,141 sites. ",
         "Rooted on Nepenthes.\n",
         "Node labels: UFBoot/SH-aLRT | site concordance factor. ",
         "Red = sCF below 40%, i.e. the sites do not support that split."),
       caption = paste0(
         "The two subgenome clades ARE the allopolyploidy. regia_B groups ",
         "with Dionaea_B; regia_A sits basal.\n",
         "That disagreement is resolved by ASTRAL and by the rate-free ",
         "four-point, both of which place regia with Dionaea in BOTH subgenomes.")) +
  theme_tree2() +
  theme(plot.title = element_text(face="bold", size=13),
        plot.subtitle = element_text(size=8.5, colour="grey30"),
        plot.caption = element_text(size=7.5, colour="grey35", hjust=0),
        legend.position = c(0.12, 0.90),
        legend.text = element_text(size=8))
p <- p + xlim(0, max(node.depth.edgelength(tr)) * 1.45)

ggsave("DR/fig/DR05_13tip_tree.pdf", p, width=10, height=7.5)
ggsave("DR/fig/DR05_13tip_tree.png", p, width=10, height=7.5, dpi=220)
cat("wrote DR/fig/DR05_13tip_tree.{pdf,png}\n")

cat("\n--- branch lengths, longest first ---\n")
bl <- tibble(tip=tr$tip.label,
             len=tr$edge.length[match(seq_along(tr$tip.label), tr$edge[,2])])
print(as.data.frame(bl %>% arrange(desc(len)) %>% mutate(len=round(len,4))),
      row.names=FALSE)
cat("\n  regia and Dionaea have the SHORTEST branches -- they are the two\n")
cat("  slowest lineages (rates 0.308 and 0.325 vs ~1.0 for the rest).\n")
cat("  That is why rate attraction had to be excluded explicitly, which the\n")
cat("  informative-gene-tree test and the four-point control both did.\n")
