#!/usr/bin/env Rscript
# ============================================================================
# PF_control -- two tree-free tests that the blocks are inherited units
#
# 1 HOMEOLOG EXCLUSIVITY
#   At a locus where a species has three copies, AAB predicts exactly two A and
#   one B. Labels were assigned to each copy INDEPENDENTLY -- nothing ties
#   copies at the same locus together -- so under noise the labels within a
#   locus are independent draws at the observed p(A), giving AAA about 30% of
#   the time. An excess of AAB is something noise cannot produce.
#
# 2 BOUNDARY CONCORDANCE ACROSS SPECIES
#   A real A/B boundary sits at a fixed point in the ancestral genome, so
#   different species should place their boundaries at the same Nepenthes
#   position. Artefactual boundaries fall wherever noise puts them and should
#   be independent between species. Null: boundary positions shuffled within
#   each species, keeping their number.
#
# Neither test uses a tree, and neither uses the block call to score the genes
# that produced it.
#
# IN  DR/out/DR02_propagated_blocks.csv, DR02_gene_delta.csv, DR/locus_meta.tsv
# OUT presentation_figures/PFcontrol.{pdf,png}, PFcontrol_stats.csv
# ============================================================================
suppressMessages(library(data.table))
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
set.seed(5)
OUT <- "presentation_figures"; dir.create(OUT, showWarnings=FALSE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }
A_COL <- "#1D9E75"; B_COL <- "#D85A30"; A_DK <- "#0F6E56"; NUL <- "#B4B2A9"; N_COL <- "#6E6D69"
NPERM <- 999

hr("0. LOAD")
P <- fread("DR/out/DR02_propagated_blocks.csv")[label %in% c("A","B") & is.finite(mb)]
gk <- if ("gene" %in% names(P)) "gene" else "tip"
G <- fread("DR/out/DR02_gene_delta.csv")[is.finite(delta)]
cat("  propagated:", format(nrow(P), big.mark=","), "| gene_delta cols:",
    paste(names(G), collapse=", "), "\n")

hr("1. HOMEOLOG EXCLUSIVITY")
if (!"locus" %in% names(G)) die("no locus column in DR02_gene_delta.csv")
L <- merge(G[, .(locus, genome, gene)],
           P[, .(genome, gene=get(gk), label)], by=c("genome","gene"))
cat(sprintf("  copies with a locus and a label: %s | loci %s\n",
    format(nrow(L), big.mark=","), format(uniqueN(L$locus), big.mark=",")))
L3 <- L[, .N, by=.(locus, genome)][N==3]
cat(sprintf("  (locus, species) with exactly three copies: %s\n",
    format(nrow(L3), big.mark=",")))
if (!nrow(L3)) die("no three-copy loci -- cannot run this test")
T3 <- merge(L, L3[, .(locus, genome)], by=c("locus","genome"))
CMB <- T3[, .(nA = sum(label=="A")), by=.(locus, genome)]
obs <- CMB[, .N, by=nA][order(nA)]
pA  <- mean(T3$label == "A")
exp <- data.table(nA=0:3, p=dbinom(0:3, 3, pA))
exp[, expected := p * nrow(CMB)]
R1 <- merge(exp[, .(nA, expected)], obs, by="nA", all.x=TRUE)
setnafill(R1, fill=0, cols="N"); setnames(R1, "N", "observed")
R1[, combo := c("BBB","ABB","AAB","AAA")[nA+1]]
R1[, ratio := observed/pmax(expected, 0.5)]
cat(sprintf("  p(A) among these copies: %.3f\n", pA))
print(R1[, .(combo, observed, expected=round(expected,1), ratio=round(ratio,2))],
      row.names=FALSE)
cs <- suppressWarnings(chisq.test(R1$observed, p=R1$expected/sum(R1$expected)))
cat(sprintf("  chi-square vs independent labels: X2 = %.1f, df = %d, p = %.3g\n",
    cs$statistic, cs$parameter, cs$p.value))
cat(sprintf("  AAB observed %d vs %.0f expected -> %.2fx\n",
    R1$observed[R1$combo=="AAB"], R1$expected[R1$combo=="AAB"],
    R1$ratio[R1$combo=="AAB"]))
cat(sprintf("  all-same (AAA + BBB) observed %d vs %.0f expected -> %.2fx\n",
    sum(R1$observed[R1$combo %in% c("AAA","BBB")]),
    sum(R1$expected[R1$combo %in% c("AAA","BBB")]),
    sum(R1$observed[R1$combo %in% c("AAA","BBB")]) /
    max(sum(R1$expected[R1$combo %in% c("AAA","BBB")]), 0.5)))

hr("2. BOUNDARY CONCORDANCE")
if (!"region" %in% names(P)) die("no region column")
bnd <- P[order(genome, chr, mb), {
  lb <- label; n <- .N
  i <- which(head(lb,-1) != tail(lb,-1))
  if (!length(i)) NULL else .(mb_b = (mb[i] + mb[i+1])/2, region = region[i]) },
  by=.(genome, chr)]
cat(sprintf("  A/B boundaries found: %d across %d chromosomes\n",
    nrow(bnd), uniqueN(bnd[, paste(genome, chr)])))
if (nrow(bnd) < 20) die("too few boundaries for this test")
print(bnd[, .N, by=genome][order(-N)], row.names=FALSE)
cat("\n  NOTE this test needs each boundary's ANCESTRAL (Nepenthes) coordinate.\n")
cat("  DR02_propagated_blocks.csv carries region but not a Nepenthes position,\n")
cat("  so the cross-species version needs combBed. Reporting per-region counts\n")
cat("  for now -- if several species place a boundary in the SAME region more\n")
cat("  often than chance, that is already suggestive.\n\n")
rb <- bnd[, .(species = uniqueN(genome), boundaries = .N), by=region][order(-species)]
print(rb, row.names=FALSE)
nsp <- uniqueN(bnd$genome)
perm <- replicate(NPERM, {
  s <- bnd[, .(region = sample(unique(P$region), .N, replace=TRUE)), by=genome]
  max(s[, .(k = uniqueN(genome)), by=region]$k) })
cat(sprintf("\n  regions where >= 3 species place a boundary: observed %d\n",
    sum(rb$species >= 3)))
cat(sprintf("  most species sharing one region: observed %d | shuffled median %d, 95th %d\n",
    max(rb$species), as.integer(median(perm)), as.integer(quantile(perm, .95))))
cat(sprintf("  p = %.3f\n", (1 + sum(perm >= max(rb$species)))/(NPERM+1)))

fwrite(R1, file.path(OUT,"PFcontrol_stats.csv"))

hr("3. FIGURE")
draw <- function(){
  par(mar=c(4.6,5.4,4.0,1.6))
  o <- R1[order(nA)]
  x <- barplot(rbind(o$observed, o$expected), beside=TRUE, border=NA,
               col=c(A_COL, NUL), names.arg=o$combo, las=1,
               cex.names=1.15, cex.axis=0.95,
               ylim=c(0, max(c(o$observed, o$expected))*1.22))
  mtext("three-copy loci", side=2, line=3.8, cex=1.0)
  mtext("how the three copies at one locus are labelled", side=1, line=2.8, cex=1.0)
  for (i in seq_len(nrow(o)))
    if (o$observed[i] > 0)
      text(x[1,i], o$observed[i], sprintf("%.1fx", o$ratio[i]),
           pos=3, cex=0.95, col=A_DK, font=2)
  legend("topleft", bty="n", cex=1.0, fill=c(A_COL, NUL), border=NA,
         legend=c("observed", "if the labels were independent"))
  mtext("Copies at the same locus are not labelled independently",
        side=3, adj=0, line=2.2, cex=1.2, font=2)
  mtext(sprintf("%s three-copy loci | nothing in the method ties copies at a locus together",
        format(sum(o$observed), big.mark=",")),
        side=3, adj=0, line=0.8, cex=0.88, col=N_COL) }
W <- 8.5; H <- 5.4
pdf(file.path(OUT,"PFcontrol.pdf"), width=W, height=H, useDingbats=FALSE); draw(); dev.off()
png(file.path(OUT,"PFcontrol.png"), width=W*160, height=H*160, res=160); draw(); dev.off()
cat("  wrote PFcontrol.pdf / .png\n")
