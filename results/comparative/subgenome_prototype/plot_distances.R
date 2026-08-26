suppressMessages({library(data.table); library(ggplot2)})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GRN <- "#1b9e77"; ORG <- "#d95f02"; PUR <- "#7570b3"
WHITE <- theme(plot.background=element_rect(fill="white", colour=NA))

save_both <- function(plot, stem, w, h) {
  pp <- normalizePath(file.path(getwd(), paste0(stem, ".pdf")), mustWork=FALSE)
  ggsave(pp, plot, width=w, height=h); cat("WROTE:", pp, "\n")
  gp <- normalizePath(file.path(getwd(), paste0(stem, ".png")), mustWork=FALSE)
  ok <- FALSE
  if (requireNamespace("ragg", quietly=TRUE))
    ok <- tryCatch({ ggsave(gp, plot, width=w, height=h, dpi=150, device=ragg::agg_png); TRUE },
                   error=function(e) FALSE)
  if (!ok) ok <- tryCatch({ ggsave(gp, plot, width=w, height=h, dpi=150,
                                   device=grDevices::png, type="cairo"); TRUE },
                          error=function(e) FALSE)
  if (ok) cat("WROTE:", gp, "\n")
}

d <- fread("distances.csv")
num <- grep("^(pdist|dN|dS|naa|alnlen)", names(d), value=TRUE)
d[, (num) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols=num]
cat("triplets:", nrow(d), " | arrangements:\n"); print(d[, .N, by=arrangement])

# ---- dS saturation diagnostic ----
# dS is Jukes-Cantor corrected. Invert it to recover the RAW fraction of
# synonymous sites that differ: p = 3/4 (1 - exp(-4d/3)). Random sequences
# differ at 3/4 of sites, so p=0.75 is the ceiling and the correction becomes
# unstable near it. Var(d) inflation = 1/(1-4p/3)^2.
sat <- rbindlist(lapply(c("dio1_dio2","nep_dio1","nep_dio2"), function(cmp) {
  v <- d[[paste0("dS_", cmp)]]; vv <- v[is.finite(v)]
  md <- median(vv); ps <- 0.75*(1-exp(-4*md/3))
  data.table(comparison=cmp, n=length(vv), median_dS=round(md,4),
             implied_pS=round(ps,4), pct_of_ceiling=round(100*ps/0.75,1),
             var_inflation=round(1/(1-4*ps/3)^2,2), n_NaN=sum(!is.finite(v)))
}))
cat("\n=== dS saturation diagnostic ===\n"); print(sat)
cat("(pct_of_ceiling <70% = dS trustworthy; var_inflation = how much JC magnifies noise)\n")
fwrite(sat, "saturation_table.csv"); cat("WROTE:", normalizePath("saturation_table.csv", mustWork=FALSE), "\n")

cat("\n=== medians ===\n")
print(d[, .(dio1_dio2=median(dS_dio1_dio2,na.rm=TRUE),
            nep_dio1 =median(dS_nep_dio1, na.rm=TRUE),
            nep_dio2 =median(dS_nep_dio2, na.rm=TRUE))])

# ---- PER-CHROMOSOME-PAIR PAIRED TEST = phasing ----
# Genome-wide, dio1/dio2 are arbitrary labels. But WITHIN a chromosome pair,
# sg1/sg2 is consistent, so a paired test per pair asks whether one homeolog is
# systematically further from Nepenthes. A consistent sign across pairs assigns
# each chromosome to a subgenome -- phasing without SubPhaser or an outgroup.
ph <- rbindlist(lapply(sort(unique(d[chrpair!="mixed", chrpair])), function(cp) {
  s <- d[chrpair==cp]
  a <- s$dS_nep_dio1; b <- s$dS_nep_dio2; ok <- is.finite(a)&is.finite(b)
  if (sum(ok) < 15) return(NULL)
  w <- suppressWarnings(wilcox.test(a[ok], b[ok], paired=TRUE))
  data.table(chrpair=cp, n=sum(ok),
             chr_dio1=s$chr1[1], chr_dio2=s$chr2[1],
             median_diff=round(median(a[ok]-b[ok]),4),
             p=signif(w$p.value,3),
             faster=ifelse(median(a[ok]-b[ok])>0, s$chr1[1], s$chr2[1]))
}))
if (nrow(ph)) {
  ph[, p_adj := signif(p.adjust(p, "BH"),3)]
  cat("\n=== per-chromosome-pair paired test: which homeolog is FURTHER from Nepenthes? ===\n")
  print(ph)
  cat("(consistent 'faster' chromosome across pairs = a subgenome-level rate difference)\n")
  fwrite(ph, "phasing_by_chrpair.csv")
  cat("WROTE:", normalizePath("phasing_by_chrpair.csv", mustWork=FALSE), "\n")
}

# ---- A: homeolog divergence ----
m1 <- melt(d[, .(pdist_dio1_dio2, dS_dio1_dio2)],
           measure.vars=c("pdist_dio1_dio2","dS_dio1_dio2"))
m1[, variable := factor(variable, labels=c("protein p-distance","dS"))]
p1 <- ggplot(m1[is.finite(value)], aes(value)) +
  geom_histogram(bins=50, fill=PUR, colour="white", linewidth=0.15) +
  facet_wrap(~variable, scales="free") +
  labs(title="A.  Divergence between the two Dionaea homeologs",
       subtitle="All pairs coalesce at the same duplication, so dS (a clock) should be TIGHT.\nProtein distance varies with per-gene selective constraint and is expected to be BROAD.",
       x=NULL, y="gene pairs") +
  theme_minimal(base_size=12) +
  theme(strip.text=element_text(face="bold"), plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(size=9, colour="grey35")) + WHITE

# ---- B: min vs max (labels arbitrary genome-wide) ----
d[, `:=`(dS_lo=pmin(dS_nep_dio1, dS_nep_dio2), dS_hi=pmax(dS_nep_dio1, dS_nep_dio2))]
p2 <- ggplot(d[is.finite(dS_lo)], aes(dS_lo, dS_hi)) +
  geom_abline(slope=1, intercept=0, colour="grey55", linetype="dashed") +
  geom_point(alpha=0.4, size=1.6, colour=ORG) + coord_equal() +
  labs(title="B.  Distance to Nepenthes: nearer vs further homeolog",
       subtitle="Labels are arbitrary genome-wide, so plotted as min vs max.\nDistance above the line = how unequal the two copies are.",
       x="dS, nearer homeolog", y="dS, further homeolog") +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(size=9, colour="grey35")) + WHITE

# ---- C: per-chromosome-pair paired differences ----
p3 <- ggplot(d[chrpair!="mixed" & is.finite(dS_nep_dio1) & is.finite(dS_nep_dio2)],
             aes(dS_nep_dio1 - dS_nep_dio2)) +
  geom_vline(xintercept=0, colour="grey40", linetype="dashed") +
  geom_histogram(bins=35, fill=GRN, colour="white", linewidth=0.15) +
  facet_wrap(~chrpair, scales="free_y", ncol=4) +
  labs(title="C.  Paired difference, per chromosome pair",
       subtitle="Within a pair the sg1/sg2 labels ARE consistent. A histogram shifted off zero means one homeolog is systematically further from Nepenthes.",
       x="dS(Nep, sg1) - dS(Nep, sg2)", y="gene pairs") +
  theme_minimal(base_size=12) +
  theme(strip.text=element_text(face="bold", size=9), plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(size=9, colour="grey35")) + WHITE

# ---- D: dS vs protein distance ----
m4 <- rbind(d[, .(pd=pdist_dio1_dio2, dS=dS_dio1_dio2, cmp="homeolog vs homeolog")],
            d[, .(pd=pdist_nep_dio1,  dS=dS_nep_dio1,  cmp="Nepenthes vs homeolog 1")],
            d[, .(pd=pdist_nep_dio2,  dS=dS_nep_dio2,  cmp="Nepenthes vs homeolog 2")])
p4 <- ggplot(m4[is.finite(pd)&is.finite(dS)], aes(pd, dS, colour=cmp)) +
  geom_point(alpha=0.35, size=1.3) +
  geom_smooth(method="loess", se=FALSE, linewidth=0.7, formula=y~x) +
  scale_colour_manual(values=c("homeolog vs homeolog"=PUR, "Nepenthes vs homeolog 1"=ORG,
                               "Nepenthes vs homeolog 2"=GRN), name=NULL) +
  labs(title="D.  dS vs protein distance, per gene pair",
       subtitle="These axes need not correlate: dS tracks TIME, protein distance tracks CONSTRAINT.\nWhat matters is that homeologs (purple) sit below both Nepenthes comparisons across the whole range.",
       x="protein p-distance", y="dS") +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom", plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(size=9, colour="grey35")) + WHITE

CAP <- sprintf(paste0("Nepenthes gracilis dominant subgenome vs Dionaea muscipula, genome-wide. %d triplets: 1 Nepenthes gene on a _dom chromosome + exactly 2 Dionaea copies in DIFFERENT tandem arrays\n",
  "(GENESPACE arrayID), so homeologs translocated onto one chromosome are kept and only true tandem duplicates dropped: %d different_chr, %d same_chr_translocated.\n",
  "MAFFT L-INS-i on peptides, back-translated with pal2nal; dN/dS by Nei-Gojobori (multi-hit codons approximated). Retained duplicates only - a fractionation-surviving, biased sample.\n",
  "NOTE: none of this separates allo- from autopolyploidy. Under a clock both homeologs are equidistant from Nepenthes regardless, so this detects RATE differences only."),
  nrow(d), sum(d$arrangement=="different_chr"), sum(d$arrangement=="same_chr_translocated"))

if (requireNamespace("patchwork", quietly=TRUE)) {
  suppressMessages(library(patchwork))
  p <- (p1 / (p2 | p4) / p3) + plot_layout(heights=c(1, 1.15, 1.3)) +
    plot_annotation(title="Dionaea homeolog divergence, relative to Nepenthes", caption=CAP,
      theme=theme(plot.title=element_text(face="bold", size=16),
                  plot.caption=element_text(size=6.8, colour="grey45", hjust=0, lineheight=1.2),
                  plot.background=element_rect(fill="white", colour=NA)))
  save_both(p, "subgenome_distances", 13, 15)
} else {
  for (i in 1:4) save_both(get(paste0("p", i)), sprintf("subgenome_%s", LETTERS[i]), 9, 6)
}
