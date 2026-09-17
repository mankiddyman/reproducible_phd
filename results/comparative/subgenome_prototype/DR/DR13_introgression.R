#!/usr/bin/env Rscript
# ============================================================================
# DR13 — IS THE B ASYMMETRY INTROGRESSION? THREE TESTS
#
# THE PROBLEM
#   DR11 (topology counting) : B D = +0.296 p = 7e-14 | A D = -0.014 p = 0.62
#   DR12 (site patterns)     : B D = 0.096 mean, 1/4 significant
#                              A D = -0.088 to -0.156, ALL FOUR NEGATIVE
#   They disagree, and DR12's A "control" is not clean, so neither can be
#   read at face value.
#
# WHY DR12 WAS WRONG
#   Its block jackknife cut the CONCATENATED matrix into 50 contiguous slices.
#   Sites within a locus share a gene tree, so blocks must be LOCI. Slices
#   straddle locus boundaries, the blocks are not independent, and the
#   standard errors are wrong.
#
# THE THREE TESTS HERE
#   1. LOCUS-LEVEL D. ABBA/BABA counted per locus, jackknifed over loci.
#      This is the standard calibration and it respects linkage.
#   2. INTERNAL CONTROLS. Quartets entirely within core Drosera, where no
#      introgression is expected. If those also deviate, the matrix has a
#      systematic problem and nothing else here is interpretable.
#   3. SPATIAL CLUSTERING. Introgressed material arrives as CONTIGUOUS
#      TRACTS, so loci with excess ABBA should cluster along the chromosome.
#      Scattered -> systematic error or ancestral structure.
#      Clustered -> introgression, and the tract extent gives the genome
#      fraction, which is the number D cannot provide.
#
# IN   DR/tree/genes/*.fna   per-locus alignments, tips named by subgenome unit
# OUT  DR/out/DR13_locus_counts.csv, DR13_dstat.csv, DR13_tracts.csv
# FIG  DR/fig/DR13_1_dstat.pdf, DR13_2_clustering.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))
set.seed(1); NPERM <- 999L
CORE <- c("binata","paradoxa","scorpioides","capensis")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

readfa <- function(p) {
  x <- readLines(p, warn=FALSE); h <- grep("^>", x)
  if (!length(h)) return(NULL)
  nm <- sub("^>", "", sub(" .*$", "", x[h]))
  st <- h+1; en <- c(h[-1]-1, length(x))
  setNames(lapply(seq_along(h), function(i)
    strsplit(toupper(paste0(x[st[i]:en[i]], collapse="")), "")[[1]]), nm)
}

hr("1. COUNT ABBA / BABA PER LOCUS")
fs <- list.files("DR/tree/genes", pattern="[.]fna$", full.names=TRUE)
cat(sprintf("  per-locus alignments: %d\n", length(fs)))
ALN <- setNames(lapply(fs, readfa), sub("[.]fna$", "", basename(fs)))
ALN <- ALN[!vapply(ALN, is.null, logical(1))]
cat(sprintf("  parsed: %d\n", length(ALN)))

counts <- function(P1, P2, P3, O) {
  bind_rows(lapply(names(ALN), function(l) {
    s <- ALN[[l]]
    if (!all(c(P1,P2,P3,O) %in% names(s))) return(NULL)
    a <- s[[P1]]; b <- s[[P2]]; c3 <- s[[P3]]; oo <- s[[O]]
    n <- min(lengths(list(a,b,c3,oo)))
    a<-a[1:n]; b<-b[1:n]; c3<-c3[1:n]; oo<-oo[1:n]
    ok <- a %in% c("A","C","G","T") & b %in% c("A","C","G","T") &
          c3 %in% c("A","C","G","T") & oo %in% c("A","C","G","T")
    if (sum(ok) < 50) return(NULL)
    a<-a[ok]; b<-b[ok]; c3<-c3[ok]; oo<-oo[ok]
    tibble(locus=l, sites=length(a),
           ABBA=sum(a==oo & b==c3 & b!=oo),
           BABA=sum(b==oo & a==c3 & a!=oo))
  }))
}
jack <- function(d) {
  A <- sum(d$ABBA); B <- sum(d$BABA)
  if (A+B < 50) return(NULL)
  D <- (A-B)/(A+B); n <- nrow(d)
  # leave-one-LOCUS-out: loci are the independent unit, not arbitrary slices
  jk <- vapply(seq_len(n), function(i) {
    ai <- A - d$ABBA[i]; bi <- B - d$BABA[i]
    if (ai+bi == 0) NA_real_ else (ai-bi)/(ai+bi)
  }, numeric(1))
  jk <- jk[!is.na(jk)]; m <- length(jk)
  se <- sqrt((m-1)/m * sum((jk - mean(jk))^2))
  tibble(loci=n, ABBA=A, BABA=B, D=D, SE=se, Z=D/se,
         p=2*pnorm(-abs(D/se)))
}

hr("2. THE FOCAL TEST — regia vs Dionaea against each core species")
FOC <- bind_rows(lapply(c("A","B"), function(sg)
  bind_rows(lapply(CORE, function(cs) {
    d <- counts(paste0("Dionaea_",sg), paste0("regia_",sg),
                paste0(cs,"_",sg), "Nepenthes")
    if (is.null(d) || !nrow(d)) return(NULL)
    j <- jack(d); if (is.null(j)) return(NULL)
    j %>% mutate(subgenome=sg, P3=cs, .before=1)
  }))))
print(as.data.frame(FOC %>% mutate(D=round(D,4), SE=round(SE,4),
  Z=round(Z,2), p=signif(p,3))), row.names=FALSE)
cat("\n  D > 0: regia shares more derived alleles with core Drosera.\n")
cat("  D < 0: Dionaea does.\n")

hr("3. INTERNAL CONTROLS — quartets with no expected introgression")
cat("  All three ingroup taxa inside core Drosera. If these deviate, the\n")
cat("  matrix has a systematic bias and section 2 is uninterpretable.\n\n")
CTL <- bind_rows(lapply(c("A","B"), function(sg)
  bind_rows(lapply(list(c("binata","paradoxa","capensis"),
                        c("binata","scorpioides","capensis"),
                        c("paradoxa","scorpioides","capensis"),
                        c("binata","paradoxa","scorpioides")), function(tr) {
    d <- counts(paste0(tr[1],"_",sg), paste0(tr[2],"_",sg),
                paste0(tr[3],"_",sg), "Nepenthes")
    if (is.null(d) || !nrow(d)) return(NULL)
    j <- jack(d); if (is.null(j)) return(NULL)
    j %>% mutate(subgenome=sg, quartet=paste(tr, collapse="/"), .before=1)
  }))))
print(as.data.frame(CTL %>% mutate(D=round(D,4), SE=round(SE,4),
  Z=round(Z,2), p=signif(p,3))), row.names=FALSE)
cat("\n  |Z| < 2 throughout -> the method is clean and section 2 stands.\n")
cat("  Controls also deviating -> systematic bias, stop reading here.\n")
write_csv(bind_rows(FOC %>% mutate(type="focal"),
                    CTL %>% mutate(type="control")), "DR/out/DR13_dstat.csv")

hr("4. SPATIAL CLUSTERING — the test that says INTROGRESSION vs noise")
cat("  Introgressed DNA arrives in contiguous tracts, so loci with excess\n")
cat("  ABBA should cluster along the chromosome. Permuting locus ORDER keeps\n")
cat("  the counts and destroys the arrangement.\n\n")
LC <- counts("Dionaea_B", "regia_B", "binata_B", "Nepenthes") %>%
  mutate(excess = ABBA - BABA,
         region = sub("-.*$", "", locus),
         ord    = suppressWarnings(as.numeric(sub("^.*-g([0-9]+).*$", "\\1", locus)))) %>%
  filter(!is.na(ord)) %>% arrange(region, ord)
cat(sprintf("  loci with positions: %d across %d regions\n",
            nrow(LC), n_distinct(LC$region)))
runstat <- function(v) {
  pos <- v > 0
  if (!any(pos)) return(0L)
  r <- rle(pos); max(r$lengths[r$values])
}
CLU <- bind_rows(lapply(unique(LC$region), function(rg) {
  d <- LC %>% filter(region==rg); if (nrow(d) < 15) return(NULL)
  obs <- runstat(d$excess)
  nul <- vapply(seq_len(NPERM), function(i) runstat(sample(d$excess)), integer(1))
  tibble(region=rg, loci=nrow(d), excess_sum=sum(d$excess),
         obs_run=obs, null_p95=quantile(nul, .95),
         p=mean(nul >= obs))
}))
if (nrow(CLU)) {
  print(as.data.frame(CLU %>% mutate(p=signif(p,3))), row.names=FALSE)
  cat(sprintf("\n  regions with a run longer than chance: %d of %d\n",
              sum(CLU$p < 0.05), nrow(CLU)))
  cat("  CLUSTERED -> introgression in tracts. The run length gives the\n")
  cat("    genome fraction, which D alone cannot.\n")
  cat("  SCATTERED -> systematic error or ancestral structure, NOT tracts.\n")
  write_csv(CLU, "DR/out/DR13_tracts.csv")
} else cat("  too few positioned loci per region\n")
write_csv(LC, "DR/out/DR13_locus_counts.csv")

p1 <- ggplot(bind_rows(FOC %>% mutate(lab=paste(subgenome, P3), type="focal"),
                       CTL %>% mutate(lab=paste(subgenome, quartet), type="control")),
             aes(reorder(lab, D), D, colour=type)) +
  geom_hline(yintercept=0, colour="grey40") +
  geom_linerange(aes(ymin=D-1.96*SE, ymax=D+1.96*SE)) +
  geom_point(size=2.4) + coord_flip() +
  scale_colour_manual(values=c(focal="#D85A30", control="#378ADD"), name=NULL) +
  labs(title="DR13 - D with a locus-level jackknife",
       subtitle="controls are quartets inside core Drosera, where D should be 0",
       x=NULL, y="D (95% CI)") +
  theme_minimal(10) + theme(legend.position="top")
suppressWarnings(ggsave("DR/fig/DR13_1_dstat.pdf", p1, width=8, height=6))
p2 <- ggplot(LC, aes(ord, excess)) +
  geom_hline(yintercept=0, colour="grey50") +
  geom_col(width=0.8, fill="#378ADD") +
  facet_wrap(~region, scales="free_x", ncol=2) +
  labs(title="DR13 - ABBA minus BABA along each ancestral region",
       subtitle="runs of positive bars = introgressed tracts | scatter = noise",
       x="gene order", y="ABBA - BABA") + theme_minimal(8)
suppressWarnings(ggsave("DR/fig/DR13_2_clustering.pdf", p2, width=9, height=10))
