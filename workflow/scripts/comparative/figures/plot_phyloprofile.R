suppressMessages({library(ggtree); library(ape); library(ggplot2); library(data.table)})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/holocentricity")
outdir <- "figures"; dir.create(outdir, showWarnings=FALSE)

SHOW_CENTROMERE_EVENTS <- TRUE   # holo-ancestral reading: 2 reversions to mono. FALSE = gene events only.

# ---- species tree: Dionaea sister to Drosera, Nepenthes outgroup ----
nwk <- paste0("(Nepenthes_gracilis,(Dionaea_muscipula,(Drosera_regia,",
              "((Drosera_binata,(Drosera_paradoxa,Drosera_roseana)),",
              "(Drosera_capensis,Drosera_aliciae,Drosera_tokaiensis,Drosera_filiformis)))));")
tr <- read.tree(text=nwk)

dat <- data.table(
  label      = c("Nepenthes_gracilis","Dionaea_muscipula","Drosera_regia","Drosera_binata",
                 "Drosera_paradoxa","Drosera_roseana","Drosera_capensis","Drosera_aliciae",
                 "Drosera_tokaiensis","Drosera_filiformis"),
  nice       = c("Nepenthes gracilis","Dionaea muscipula","D. regia","D. binata",
                 "D. paradoxa","D. roseana","D. capensis","D. aliciae",
                 "D. tokaiensis","D. filiformis"),
  centromere = c("monocentric","monocentric","holocentric","monocentric","holocentric",
                 "holocentric","monocentric","monocentric","monocentric","monocentric"),
  knl2       = c("present","absent","absent","absent","absent",
                 "absent","absent","absent","absent","absent"),
  syn4       = c(2, 1, 1, 4, 2, 2, 4, 6, 6, 3),
  tpx2       = c("present","present","present","silenced","present",
                 "present","absent","absent","absent","absent")
)
dat[, syn4_lab := fifelse(label=="Nepenthes_gracilis", "2*", as.character(syn4))]

cent_cols  <- c(holocentric="#1b9e77", monocentric="#d95f02")
state_fill <- c(present="grey15", absent="white", silenced="#d95f02")
state_shape<- c(present=21, absent=21, silenced=8)

p <- ggtree(tr, size=0.9, branch.length="none") %<+% dat +
  geom_tippoint(aes(color=centromere), size=3.6) +
  geom_tiplab(aes(label=nice, color=centromere), fontface="italic", size=4.3, offset=0.15)

pd   <- as.data.table(p$data)
tipd <- merge(pd[isTip==TRUE, .(label, x, y)], dat, by="label")
ntip <- length(tr$tip.label); xmax <- max(pd$x)
BW <- 0.17
c1 <- xmax + 5.4                 # KNL2
c2 <- c1 + 1.3                   # SYN4-like bars
c3 <- c2 + 6*BW + 1.2            # TPXL6
xhi <- c3 + 1.2
tipd[, `:=`(x1=c1, x2=c2, x3=c3)]

p <- p +
  geom_point(data=tipd, aes(x=x1, y=y, shape=knl2, fill=knl2),
             size=4.6, color="grey25", stroke=1, inherit.aes=FALSE) +
  geom_segment(data=tipd[syn4>0], aes(x=x2, xend=x2+syn4*BW, y=y, yend=y, color=centromere),
               linewidth=6, inherit.aes=FALSE) +
  geom_text(data=tipd, aes(x=x2+syn4*BW+0.10, y=y, label=syn4_lab),
            size=3.5, hjust=0, color="grey25", inherit.aes=FALSE) +
  geom_point(data=tipd, aes(x=x3, y=y, shape=tpx2, fill=tpx2),
             size=4.6, color="grey25", stroke=1, inherit.aes=FALSE) +
  scale_color_manual(values=cent_cols, name="centromere") +
  scale_fill_manual(values=state_fill, name="gene state") +
  scale_shape_manual(values=state_shape, name="gene state")

# ---- branch events ----
n_drose <- getMRCA(tr, c("Dionaea_muscipula","Drosera_regia"))
n_dros  <- getMRCA(tr, c("Drosera_regia","Drosera_capensis"))
n_cap   <- getMRCA(tr, c("Drosera_capensis","Drosera_aliciae"))
n_bin   <- which(tr$tip.label=="Drosera_binata")
midx <- function(nd){ xn <- pd[node==nd, x]; xp <- pd[node==pd[node==nd, parent], x]; (xn+xp)/2 }
midy <- function(nd) pd[node==nd, y]

ev <- data.table(node=c(n_drose, n_bin, n_cap),
                 lab =c("KNL2 lost\nSYN4-like origin", "TPXL6 silenced?", "TPXL6 deleted"))
if (SHOW_CENTROMERE_EVENTS)
  ev <- rbind(ev, data.table(node=c(n_dros, n_bin, n_cap),
                             lab=c("holocentric", "→ monocentric", "→ monocentric")))
ev[, `:=`(ex=sapply(node, midx), ey=sapply(node, midy))]
ev[, ey2 := ey + c(0.30,-0.34,-0.34, 0.30, 0.30, 0.30)[seq_len(.N)]]

p <- p +
  geom_point(data=ev, aes(x=ex, y=ey), shape=23, size=3.4, fill="grey35",
             color="black", inherit.aes=FALSE) +
  geom_text(data=ev, aes(x=ex-0.12, y=ey2, label=lab), hjust=1, size=3.1,
            lineheight=0.9, color="grey20", inherit.aes=FALSE)

p <- p +
  annotate("text", x=c(c1, c2+3*BW, c3), y=ntip+1.1, size=4, fontface="bold", lineheight=0.9,
           label=c("α/βKNL2", "\u03b1-kleisin\nSYN4 clade", "TPXL6-short\n(AT5G37478)")) +
  expand_limits(x=c(0, xhi), y=c(0.2, ntip+2.1)) +
  labs(title="Centromere-machinery changes across the Droseraceae",
       caption=paste0("*Nepenthes 2 = canonical SYN4; all others = the Droseraceae-specific SYN4-like clade. Drosera = raw OG0001569 copies; Nepenthes/Dionaea counted from the kleisin tree (not in the 8-Drosera OrthoFinder run).  ",
                      "KNL2 = tblastn+SANTA(PF09133) vs genomes.  SYN4-like = raw OG0001569 copies.\n",
                      "TPXL6-short: presence/absence genome-verified; binata relic = frameshifted, stop-riddled, no ORF.  ",
                      "n=2 independent TPXL6 losses — no statistical support possible.")) +
  theme_tree(legend.position=c(0.09, 0.72)) +
  theme(plot.title=element_text(face="bold", size=17),
        plot.caption=element_text(size=7.5, color="grey40", hjust=0),
        legend.text=element_text(size=10), legend.title=element_text(size=10, face="bold"))

ggsave(file.path(outdir,"phylo_profile.pdf"), p, width=13, height=7.5)
ggsave(file.path(outdir,"phylo_profile.jpg"), p, width=13, height=7.5, dpi=200)
cat("wrote phylo_profile.{pdf,jpg}\n")
print(tipd[, .(nice, centromere, knl2, syn4, tpx2)])
