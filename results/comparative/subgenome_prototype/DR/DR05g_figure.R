#!/usr/bin/env Rscript
# ============================================================================
# DR05g — THE FIGURE
#
# PANEL A/B: ASTRAL, subgenomes A and B, rooted on Nepenthes.
#   All four ASTRAL runs (A and B, full and informative-only) give the SAME
#   rooted topology. Two independent partitions, same answer.
#
# WHY ASTRAL AND NOT CONCATENATION as the primary tree:
#   Not a general preference. At THIS node concatenation is measurably wrong.
#   A_full node 8 -- the branch separating regia_A -- has sCF 32.66 with
#   sDF2 41.61: the ALTERNATIVE resolution has more site support than the
#   topology in the tree, while bootstrap reports 100/100. Panel C shows it.
#
# BRANCH LENGTHS: ASTRAL internal branches are in COALESCENT UNITS; terminal
#   branches are not estimated. So the tree is drawn as a cladogram and the
#   key internal branch is labelled with its CU value. 0.22-0.25 CU is SHORT
#   and the reader should see that rather than only the local PP of 1.0.
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(ape); library(ggtree); library(ggplot2); library(patchwork); library(dplyr)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))

pretty_tips <- function(x) {
  x <- sub("_A$", "", sub("_B$", "", x))
  ifelse(x=="Nepenthes", "Nepenthes gracilis",
  ifelse(x=="Dionaea",   "Dionaea muscipula",
         paste0("Drosera ", x)))
}

astral_panel <- function(f, sub, title) {
  t <- read.tree(f)
  t <- root(t, outgroup="Nepenthes", resolve.root=TRUE)
  t$node.label[is.na(t$node.label) | t$node.label==""] <- ""
  lab <- pretty_tips(t$tip.label)
  key <- which(t$tip.label %in% c(paste0("regia_", sub), paste0("Dionaea_", sub)))
  mrca <- getMRCA(t, key)
  cu <- t$edge.length[which(t$edge[,2]==mrca)]
  d <- data.frame(label=t$tip.label, pretty=lab,
                  grp=ifelse(t$tip.label %in% t$tip.label[key], "snap-trap clade",
                      ifelse(t$tip.label=="Nepenthes", "outgroup", "core Drosera")))
  p <- ggtree(t, branch.length="none", size=0.7) %<+% d +
    geom_tiplab(aes(label=pretty, colour=grp), fontface=3, size=3.4, offset=0.15) +
    geom_nodepoint(size=1.6, colour="grey35") +
    scale_colour_manual(values=c("snap-trap clade"="#D85A30",
                                 "core Drosera"="#1D9E75",
                                 "outgroup"="grey45"), guide="none") +
    geom_hilight(node=mrca, fill="#D85A30", alpha=0.10, extend=2.6) +
    geom_cladelab(node=mrca, label=sprintf("%.2f CU", cu),
                  offset=-0.55, offset.text=0.06, barsize=0,
                  fontsize=2.9, textcolour="#993C1D") +
    labs(title=title,
         subtitle="local PP = 1.0 at every node") +
    xlim(0, 7) + theme_tree() +
    theme(plot.title=element_text(face="bold", size=11),
          plot.subtitle=element_text(size=8, colour="grey35"))
  list(plot=p, cu=cu)
}

A <- astral_panel("DR/tree/astral/astral_A_informative.tre", "A",
                  "A subgenome")
B <- astral_panel("DR/tree/astral/astral_B_informative.tre", "B",
                  "B subgenome")
cat(sprintf("(regia,Dionaea) branch: A = %.3f CU | B = %.3f CU\n", A$cu, B$cu))

# ---- panel C: the concatenation conflict, with sCF ----
cf <- read.table("DR/tree/concat/A_full_scf.cf.stat", header=TRUE, comment.char="#")
tc <- read.tree("DR/tree/concat/A_full.treefile")
tc <- root(tc, outgroup="Nepenthes", resolve.root=TRUE)
lab <- pretty_tips(tc$tip.label)
dc <- data.frame(label=tc$tip.label, pretty=lab)
sc <- cf %>% transmute(node=ID, scf=sCF,
                       txt=sprintf("sCF %.0f%%", sCF),
                       weak=sCF < 40)
pC <- ggtree(tc, branch.length="none", size=0.7) %<+% dc +
  geom_tiplab(aes(label=pretty), fontface=3, size=3.4, offset=0.15) +
  geom_nodepoint(size=1.6, colour="grey35")
pC <- pC %<+% sc +
  geom_text2(aes(subset=!isTip, label=txt, colour=weak),
             hjust=1.15, vjust=-0.7, size=2.7) +
  scale_colour_manual(values=c("FALSE"="grey35","TRUE"="#D85A30"), guide="none") +
  labs(title="Concatenated tree, A subgenome",
       subtitle="all nodes UFBoot 100/100 | sCF shows the sites disagree at the red node") +
  xlim(0, 7) + theme_tree() +
  theme(plot.title=element_text(face="bold", size=11),
        plot.subtitle=element_text(size=8, colour="grey35"))

fig <- (A$plot | B$plot) / pC +
  plot_annotation(
    title = "Drosera regia falls outside core Drosera, with Dionaea",
    subtitle = sprintf(paste0("ASTRAL on %s informative gene trees per subgenome. ",
      "Both subgenomes give the same rooted topology.\n",
      "Rate-free four-point: regia+Dionaea in 52.3%% of loci vs 6.9-11.4%% for ",
      "any other species in the same slot."), "500+"),
    theme = theme(plot.title=element_text(face="bold", size=13),
                  plot.subtitle=element_text(size=9, colour="grey30")))
ggsave("DR/fig/DR05_MAIN_tree.pdf", fig, width=11, height=9)
ggsave("DR/fig/DR05_MAIN_tree.png", fig, width=11, height=9, dpi=200)
cat("wrote DR/fig/DR05_MAIN_tree.{pdf,png}\n")

# ---- supplementary: the rate-free four-point control ----
ctl <- read.csv("DR/out/DR05f_control.csv")
ps <- ggplot(ctl, aes(reorder(focal, focal_Dionaea), focal_Dionaea,
                      fill=focal=="regia")) +
  geom_hline(yintercept=1/3, linetype="dashed", colour="grey40") +
  geom_col(width=0.62) +
  geom_text(aes(label=sprintf("%.3f", focal_Dionaea)), hjust=-0.25, size=3.2) +
  annotate("text", x=1.4, y=0.345, hjust=0, size=2.9, colour="grey35",
           label="1/3 = no signal") +
  scale_fill_manual(values=c("FALSE"="grey70","TRUE"="#D85A30"), guide="none") +
  coord_flip(ylim=c(0,0.62)) +
  labs(title="Rate-free four-point control",
       subtitle=paste0("Each species placed in the focal slot of ",
                       "(focal_X, Dionaea_X, other_X, Nepenthes).\n",
                       "Buneman's condition: every tip enters each of the three ",
                       "distance sums once, so lineage rate cancels."),
       x=NULL, y="fraction of loci grouping focal with Dionaea") +
  theme_minimal(11) +
  theme(plot.title=element_text(face="bold"),
        plot.subtitle=element_text(size=8, colour="grey30"))
ggsave("DR/fig/DR05_SUPP_fourpoint.pdf", ps, width=8, height=4.5)
cat("wrote DR/fig/DR05_SUPP_fourpoint.pdf\n")
