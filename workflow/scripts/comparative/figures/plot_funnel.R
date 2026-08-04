suppressMessages({library(ggplot2); library(data.table)})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/holocentricity")
outdir <- "figures"

GRN <- "#1b9e77"; ORG <- "#d95f02"; GRY <- "grey55"

d <- data.table(
  step = 1:6,
  n    = c(27693, 118, 7, 3, 1, 1),
  stage = c("orthogroups tested",
            "strict complete P/A",
            "chromosome-related",
            "coherent GO signature",
            "survives reciprocal BLAST",
            "absence verified in genomes"),
  how = c("OrthoFinder, 8 Drosera",
          "present in ALL of one centromere group,\nabsent in ALL of the other",
          "broad keyword net on GO description\n(centromere|kinetochore|cohesin|CENP|meiosis|...)",
          "GO term list internally consistent\nwith a real chromosome-biology gene",
          "no full-length homolog in the 'absent' group,\ncalibrated against each OG's own present-group range",
          "tblastn vs genomes, not annotations:\nonly the shared domain of a longer paralog remains"),
  lost = c(NA,
           "27,575 not group-specific",
           "111 no chromosome term",
           "4 FANTASIA misannotation\n(animal immune terms, stress TFs)",
           "2 orthogroup-splitting artifacts\n(homolog was there all along)",
           NA)
)
d[, ymax := -step]

# funnel width on log scale
d[, w := log10(n + 1)]
d[, w := w / max(w)]
BW <- 0.42
d[, `:=`(xmin = -w*BW, xmax = w*BW)]
d[, wide := (xmax - xmin) > 0.18]   # can the box hold its own label?

lab <- data.table(
  step = c(2,3,4,5,6),
  txt  = c("118", "7", "3", "1", "1")
)

p <- ggplot(d, aes(y=ymax)) +
  geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymax-0.32, ymax=ymax+0.32,
                fill=factor(fifelse(step==6, "keep", fifelse(step==5,"keep","mid")))),
            color="grey25", linewidth=0.4) +
  # count: inside the box when it fits, otherwise immediately right of it
  geom_text(data=d[wide==TRUE], aes(x=0, label=format(n, big.mark=",")),
            size=5, fontface="bold", color="white") +
  geom_text(data=d[wide==FALSE], aes(x=xmax+0.04, label=format(n, big.mark=",")),
            size=5, fontface="bold", color="grey15", hjust=0) +
  # stage name, left
  geom_text(aes(x=-BW-0.06, label=stage), hjust=1, size=4.2, fontface="bold", color="grey15") +
  # method, right
  geom_text(aes(x=BW+0.06, label=how), hjust=0, size=3.1, color="grey40", lineheight=0.92) +
  # what was lost, between rows
  geom_text(data=d[!is.na(lost)], aes(x=-BW-0.06, y=ymax+0.52, label=paste0("\u2717  ", lost)),
            hjust=1, size=2.9, color=ORG, fontface="italic", lineheight=0.92) +
  scale_fill_manual(values=c(mid=GRY, keep=GRN), guide="none") +
  scale_x_continuous(limits=c(-1.55, 1.35), expand=c(0,0)) +
  scale_y_continuous(limits=c(-6.75, -0.35), expand=c(0,0)) +
  labs(title="One finding survived: attrition of the presence/absence screen",
       subtitle="Each row is a filter. The screen's real output is a rigorous negative \u2014 six of seven chromosome-related\ncandidates were OrthoFinder orthogroup-splitting artifacts, invisible until tested against the genomes.",
       caption=paste0("Fisher p floor at 3 holo vs 5 mono = 1/C(8,3) = 0.0179; q_BH = 1 for every candidate. ",
                      "The screen cannot produce significance \u2014 these are candidates, not biomarkers.\n",
                      "Survivor = OG0018771 = TPXL6-short (closest Arabidopsis relative AT5G37478, 178 aa). Deleted in the capensis clade; intact but unexpressed in D. binata. ",
                      "Two independent events, both on monocentric branches. Function unknown \u2014 the spindle GO is family-level transfer to a 140 aa domain-only protein.")) +
  theme_void(base_size=13) +
  theme(plot.background=element_rect(fill="white", color=NA),
        plot.title=element_text(face="bold", size=17, hjust=0, margin=margin(b=4)),
        plot.subtitle=element_text(size=10.5, color="grey35", hjust=0, lineheight=1.05, margin=margin(b=14)),
        plot.caption=element_text(size=7.5, color="grey45", hjust=0, lineheight=1.1, margin=margin(t=10)),
        plot.margin=margin(14,14,10,14))

ggsave(file.path(outdir,"attrition_funnel.pdf"), p, width=13, height=8)
ggsave(file.path(outdir,"attrition_funnel.jpg"), p, width=13, height=8, dpi=200)
cat("wrote attrition_funnel.{pdf,jpg}\n")
print(d[, .(step, n, stage)])
