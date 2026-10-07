#!/usr/bin/env Rscript
# ============================================================================
# PF_depths -- every divergence depth, in dS, computed fresh and label-free
#
# WHAT IS LABEL-FREE (no Delta, no A/B call anywhere):
#   A/B split        median dS between Dionaea's two homeologous copies. The
#                    pairing is from synteny + fractionation only.
#   Drosera-Dionaea  median dS between a Drosera copy and a Dionaea copy,
#                    pooling both Dionaea sides (so it mixes the two depths --
#                    the MINIMUM per locus is closer to the true speciation).
#   Drosera-Drosera  median dS between copies of two Drosera species.
#   vs Nepenthes     median dS to the outgroup.
#
# Every figure below is a MEDIAN with a bootstrap 95% CI. dS is substitutions
# per synonymous site, not time: a rate difference between lineages shows up
# here as a depth difference.
#
# OUT presentation_figures/PFdepths.csv
# ============================================================================
suppressMessages(library(data.table))
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
set.seed(4)
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }
ci <- function(v, B=1000) {
  if (length(v) < 20) return(c(NA,NA))
  quantile(replicate(B, median(sample(v, replace=TRUE))), c(.025,.975)) }
row <- function(lab, v) data.table(comparison=lab, n=length(v),
  median=median(v), lo=ci(v)[1], hi=ci(v)[2],
  q25=quantile(v,.25), q75=quantile(v,.75))

hr("0. INPUT")
K <- fread("DR/out/pairwise_ks.csv")[is.finite(dS) & dS >= 0 & dS < 5 & codons >= 100]
L <- fread("DR/locus_meta.tsv")
P <- fread("fractionation_by_chrpair.csv")
s1 <- P$retained_more == P$chrA
SIDE <- c(setNames(ifelse(s1,"A","B"), P$chrA), setNames(ifelse(s1,"B","A"), P$chrB))
cat(sprintf("  dS rows kept (0 <= dS < 5, codons >= 100): %s of %s\n",
    format(nrow(K), big.mark=","),
    format(nrow(fread("DR/out/pairwise_ks.csv")), big.mark=",")))
DROS <- sort(unique(L$genome[grepl("^Drosera_", L$genome)]))
CORE <- setdiff(DROS, "Drosera_regia")
cat("  core Drosera:", paste(sub("Drosera_","", CORE), collapse=", "), "\n")

R <- list()

hr("1. THE A/B SPLIT -- label-free")
dio <- L[genome=="Dionaea_muscipula"]; dio[, side := unname(SIDE[chr])]
dio <- dio[!is.na(side)]
KD <- K[sp1=="Dionaea_muscipula" & sp2=="Dionaea_muscipula"]
a <- dio$side[match(KD$seq1, dio$tip)]; b <- dio$side[match(KD$seq2, dio$tip)]
AB <- KD$dS[!is.na(a) & !is.na(b) & a != b]
WI <- KD$dS[!is.na(a) & !is.na(b) & a == b]
if (length(AB) < 50) die("too few Dionaea A-vs-B pairs")
R[[length(R)+1]] <- row("Dionaea A vs B  (= the A/B split)", AB)
if (length(WI) >= 20) R[[length(R)+1]] <- row("Dionaea within one side (control)", WI)
cat(sprintf("  A vs B: n=%d median %.4f | within a side: n=%d median %.4f\n",
    length(AB), median(AB), length(WI), if (length(WI)) median(WI) else NA))
cat("  ^ within-side should be MUCH smaller; if not, the sg1/sg2 assignment is wrong\n")

hr("2. DROSERA vs DIONAEA")
KS <- K[(sp1 %in% DROS & sp2=="Dionaea_muscipula") | (sp2 %in% DROS & sp1=="Dionaea_muscipula")]
KS[, dsp := fifelse(sp1 %in% DROS, sp1, sp2)]
for (g in c("Drosera_regia", CORE)) {
  v <- KS[dsp==g, dS]
  if (length(v) >= 20) R[[length(R)+1]] <- row(sprintf("%s vs Dionaea", sub("Drosera_","D. ",g)), v) }
R[[length(R)+1]] <- row("core Drosera vs Dionaea (pooled)", KS[dsp %in% CORE, dS])
cat("  NOTE this pools both Dionaea sides, so it sits BETWEEN the true speciation\n")
cat("  depth and the A/B depth. The per-locus MINIMUM is the better proxy:\n")
mn <- KS[, .(m = min(dS)), by=.(anchor, dsp)]
for (g in c("Drosera_regia")) {
  v <- mn[dsp==g, m]; if (length(v) >= 20)
    R[[length(R)+1]] <- row("D. regia vs Dionaea (per-locus minimum)", v) }
R[[length(R)+1]] <- row("core Drosera vs Dionaea (per-locus minimum)", mn[dsp %in% CORE, m])

hr("3. WITHIN DROSERA")
KK <- K[sp1 %in% DROS & sp2 %in% DROS & sp1 != sp2]
R[[length(R)+1]] <- row("core Drosera vs core Drosera", KK[sp1 %in% CORE & sp2 %in% CORE, dS])
R[[length(R)+1]] <- row("D. regia vs core Drosera",
                        KK[xor(sp1=="Drosera_regia", sp2=="Drosera_regia"), dS])

hr("4. TO THE OUTGROUP")
KN <- K[sp1=="Nepenthes_gracilis" | sp2=="Nepenthes_gracilis"]
KN[, other := fifelse(sp1=="Nepenthes_gracilis", sp2, sp1)]
R[[length(R)+1]] <- row("Dionaea vs Nepenthes", KN[other=="Dionaea_muscipula", dS])
R[[length(R)+1]] <- row("Drosera vs Nepenthes", KN[other %in% DROS, dS])

hr("4b. SAME-SIDE DEPTHS -- the actual speciation events")
cat("  These need A/B labels, so they are NOT label-free. Reported as depths,\n")
cat("  not as predictions. Each compares like with like: A to A, B to B.\n\n")
GD <- fread("DR/out/DR02_gene_delta.csv")[is.finite(delta)]
PB <- fread("DR/out/DR02_propagated_blocks.csv")[label %in% c("A","B")]
gk <- if ("gene" %in% names(PB)) "gene" else "tip"
LAB <- unique(merge(GD[, .(genome, gene, tip)],
                    PB[, .(genome, gene=get(gk), label)], by=c("genome","gene"))[, .(tip, label)])
SIDEMAP <- c(setNames(LAB$label, LAB$tip),
             setNames(unname(SIDE[dio$chr]), dio$tip))
K2 <- copy(K)
K2[, s_1 := unname(SIDEMAP[seq1])]; K2[, s_2 := unname(SIDEMAP[seq2])]
SS <- K2[!is.na(s_1) & !is.na(s_2) & s_1 == s_2]
cat(sprintf("  comparisons with a side on both tips: %s | same-side: %s\n",
    format(sum(!is.na(K2$s_1) & !is.na(K2$s_2)), big.mark=","),
    format(nrow(SS), big.mark=",")))
grab <- function(g1, g2, lab) {
  v <- SS[(sp1 %in% g1 & sp2 %in% g2) | (sp2 %in% g1 & sp1 %in% g2), dS]
  if (length(v) >= 20) { cat(sprintf("  %-46s n=%6d  dS %.3f\n", lab, length(v), median(v)))
                         R[[length(R)+1]] <<- row(lab, v) }
  else cat(sprintf("  %-46s too few (%d)\n", lab, length(v))) }
grab("Drosera_regia", "Dionaea_muscipula", "regia-Dionaea speciation (same side)")
grab(CORE, "Dionaea_muscipula",            "coreDrosera-Dionaea speciation (same side)")
grab("Drosera_regia", CORE,                "regia-coreDrosera speciation (same side)")
grab(CORE, CORE,                           "within core Drosera (same side)")
cat("\n  sanity: each should be SMALLER than the A/B split (0.624) if that split\n")
cat("  predates every speciation. Rate differences can still break this.\n")

hr("5. ALL DEPTHS, dS")
D <- rbindlist(R)
TAB <- D$median[D$comparison=="Dionaea A vs B  (= the A/B split)"]
D[, rel_to_AB := median/TAB]
print(D[, .(comparison, n, dS=round(median,3),
            CI=sprintf("[%.3f, %.3f]", lo, hi),
            IQR=sprintf("%.2f-%.2f", q25, q75),
            vs_AB=round(rel_to_AB,2))], row.names=FALSE)
fwrite(D, file.path(OUT,"PFdepths.csv"))
cat(sprintf("\n  A/B split = %.3f dS. Everything else is expressed relative to it.\n", TAB))
cat("  A ratio below 1 means that split is MORE RECENT than the A/B split.\n")
cat("  Delta for an A copy = (its own-side depth / the A/B depth) - 1, so the\n")
cat("  further below 1 a species sits, the further its Delta from zero.\n")
