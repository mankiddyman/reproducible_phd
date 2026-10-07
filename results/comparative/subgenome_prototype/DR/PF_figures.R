#!/usr/bin/env Rscript
# ============================================================================
# PF -- presentation figures, all from real data
#   PF1  anchor attrition cascade      (what fraction of the proteome is usable)
#   PF2  fractionation: real retention strip + observed-vs-null + per-pair scatter
#   PF3  shuffle control: real A/B tract order vs label-shuffled null
# Run in the smk env (ggplot2 + data.table), NOT the genespace module.
# OUT presentation_figures/PF{1,2,3}_*.pdf/.png
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=200)
set.seed(1)
GS   <- file.path(dirname(getwd()), "genespace")
OUT  <- "presentation_figures"
dir.create(OUT, showWarnings=FALSE)
A_COL <- "#1D9E75"; B_COL <- "#D85A30"; NULL_COL <- "#888780"
hr  <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
die <- function(...) { cat("\n  *** ", sprintf(...), " ***\n", sep=""); quit(status=1) }
th <- theme_minimal(12) + theme(
  panel.grid.minor=element_blank(), panel.grid.major=element_line(linewidth=0.2),
  plot.title=element_text(face="bold", size=14), plot.subtitle=element_text(size=10, colour="grey30"),
  axis.title=element_text(size=11), legend.position="bottom", legend.title=element_blank())

# --- 0. GATE: inspect every input before using it ---------------------------
hr("0. INPUTS")
CB <- file.path(GS,"results","combBed.txt"); PAIRF <- "fractionation_by_chrpair.csv"
PROP <- "DR/out/DR02_propagated_blocks.csv"
for (f in c(CB, PAIRF, PROP)) if (!file.exists(f)) die("missing input: %s", f)

bed <- fread(CB)
cat("\n-- combBed:", nrow(bed), "rows\n"); print(names(bed))
pairs <- fread(PAIRF)
cat("\n-- fractionation_by_chrpair.csv:", nrow(pairs), "rows\n"); print(names(pairs)); print(pairs)
prop <- fread(PROP)
cat("\n-- DR02_propagated_blocks.csv:", nrow(prop), "rows\n"); print(names(prop)); print(head(prop,3))

need_bed <- c("genome","chr","id","ofID","ord","start","end")
miss <- setdiff(need_bed, names(bed)); if (length(miss)) die("combBed lacks: %s", paste(miss, collapse=", "))
ogcol <- intersect(c("globOG","og","globHOG"), names(bed))[1]
if (is.na(ogcol)) die("combBed has no orthogroup column")
cat("\n  orthogroup column in use:", ogcol, "\n")
if (!"isArrayRep" %in% names(bed)) cat("  NOTE: no isArrayRep -- tandem filter will be skipped\n")
if (!"noAnchor"  %in% names(bed)) cat("  NOTE: no noAnchor  -- syntenic filter will be skipped\n")

pcol <- function(pat) grep(pat, names(pairs), value=TRUE, ignore.case=TRUE)[1]
cA <- pcol("^chrA$"); cB <- pcol("^chrB$")
rM <- pcol("retained_more"); rL <- pcol("retained_less|retained_fewer")
cat("  pair columns detected -- chrA:", cA, "| chrB:", cB, "| more:", rM, "| less:", rL, "\n")
if (any(is.na(c(cA,cB,rM)))) die("fractionation_by_chrpair.csv: expected chrA/chrB/retained_more; see names above")
if (is.na(rL)) cat("  NOTE: no retained_less column -- scatter will derive it\n")

pcol2 <- function(pat) grep(pat, names(prop), value=TRUE, ignore.case=TRUE)[1]
gcol <- pcol2("^gene$|^id$|^ofID$"); lcol <- pcol2("^label$"); rcol2 <- pcol2("^region$")
mcol <- pcol2("^mid$|^mb$|^start$|^pos$")
cat("  propagated columns -- gene:", gcol, "| label:", lcol, "| region:", rcol2, "| pos:", mcol, "\n")
if (any(is.na(c(gcol,lcol,rcol2,mcol)))) die("DR02_propagated_blocks.csv: missing expected columns; see names above")

# A = the more-retained Dionaea subgenome, by definition
s1 <- pairs[[rM]] == pairs[[cA]]
DIO <- c(setNames(fifelse(s1,"A","B"), pairs[[cA]]), setNames(fifelse(s1,"B","A"), pairs[[cB]]))
cat("\n  Dionaea chromosomes assigned:", length(DIO), "| A:", sum(DIO=="A"), "B:", sum(DIO=="B"), "\n")

# ============================ PF1: ATTRITION ================================
hr("PF1. ANCHOR ATTRITION")
nep <- bed[genome == "Nepenthes_gracilis"]
dio <- bed[genome == "Dionaea_muscipula"]
if (!nrow(nep) || !nrow(dio)) die("no Nepenthes (%d) or Dionaea (%d) rows in combBed", nrow(nep), nrow(dio))
dio[, sg := DIO[chr]]
cat("  Dionaea genes with a subgenome side:", sum(!is.na(dio$sg)), "of", nrow(dio), "\n")

dsum <- dio[!is.na(sg), .(nA = sum(sg=="A"), nB = sum(sg=="B")), by = c(ogcol)]
nep[, og_ := get(ogcol)]
steps <- list(); keep <- rep(TRUE, nrow(nep))
steps[["all annotated genes"]] <- sum(keep)
keep <- keep & nep$og_ %in% dsum[[ogcol]]
steps[["in an orthogroup with Dionaea"]] <- sum(keep)
ok11 <- dsum[nA == 1 & nB == 1][[ogcol]]
keep <- keep & nep$og_ %in% ok11
steps[["one copy per Dionaea subgenome"]] <- sum(keep)
if ("isArrayRep" %in% names(nep)) {
  keep <- keep & !(nep$isArrayRep %in% FALSE)
  steps[["not a tandem-array duplicate"]] <- sum(keep) }
if ("noAnchor" %in% names(nep)) {
  keep <- keep & !(nep$noAnchor %in% TRUE)
  steps[["inside a syntenic block"]] <- sum(keep) }
A1 <- data.table(step = factor(names(steps), levels = rev(names(steps))),
                 n = unlist(steps))
A1[, pct := 100*n/max(n)]
print(A1, row.names=FALSE)
ANCH <- nep$ofID[keep]
cat(sprintf("\n  final anchors: %d (%.1f%% of the Nepenthes proteome)\n", length(ANCH), 100*length(ANCH)/nrow(nep)))
if (!length(ANCH)) die("no anchors survived -- check the filters above")

p1 <- ggplot(A1, aes(n, step)) +
  geom_col(aes(fill = pct), width = 0.62) +
  geom_text(aes(label = sprintf("%s  (%.0f%%)", format(n, big.mark=","), pct)),
            hjust = -0.08, size = 3.6, colour = "grey20") +
  scale_fill_gradient(low = "#9FE1CB", high = "#0F6E56", guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(title = "What fraction of the proteome can serve as an anchor?",
       subtitle = sprintf("Nepenthes genes surviving each filter. Everything downstream rests on the final %s.",
                          format(length(ANCH), big.mark=",")),
       x = "genes", y = NULL) + th + theme(legend.position="none")
ggsave(file.path(OUT,"PF1_anchor_attrition.pdf"), p1, width=9, height=4.2)
ggsave(file.path(OUT,"PF1_anchor_attrition.png"), p1, width=9, height=4.2, dpi=200)
cat("  wrote PF1\n")

# ======================= PF2: FRACTIONATION =================================
hr("PF2. FRACTIONATION")
# per anchor: is the OG present on the A-side and on the B-side of each pair?
dio2 <- dio[!is.na(sg)]
nepA <- nep[keep, .(ofID, chr, ord, og_)]
setorder(nepA, chr, ord)
prs <- data.table(chrA = pairs[[cA]], chrB = pairs[[cB]])
prs[, `:=`(sgA = DIO[chrA], sgB = DIO[chrB])]

retain <- rbindlist(lapply(seq_len(nrow(prs)), function(i) {
  oa <- dio2[chr == prs$chrA[i], unique(get(ogcol))]
  ob <- dio2[chr == prs$chrB[i], unique(get(ogcol))]
  sub <- nepA[og_ %in% union(oa, ob)]
  if (!nrow(sub)) return(NULL)
  data.table(pair = i, chrA = prs$chrA[i], chrB = prs$chrB[i],
             ofID = sub$ofID, nep_chr = sub$chr, nep_ord = sub$ord,
             onA = sub$og_ %in% oa, onB = sub$og_ %in% ob) }))
if (!nrow(retain)) die("no anchors mapped to any chromosome pair")
ps <- retain[, .(anchors = .N, kept_A = sum(onA), kept_B = sum(onB)), by = .(pair, chrA, chrB)]
ps[, `:=`(more = pmax(kept_A, kept_B), less = pmin(kept_A, kept_B))]
ps[, asym := (more - less)/(more + less)]
print(ps[order(-asym)], row.names=FALSE)
cat(sprintf("\n  pairs where A (more-retained side) keeps more: %d of %d\n",
            sum(ps$kept_A >= ps$kept_B), nrow(ps)))
bt <- binom.test(sum(ps$kept_A > ps$kept_B), nrow(ps), 0.5)
cat(sprintf("  sign test: p = %.3g | median asymmetry %.3f\n", bt$p.value, median(ps$asym)))

# the strip: the pair with the most anchors, genes in Nepenthes order
tp <- ps[which.max(anchors)]
cat(sprintf("\n  strip drawn from pair %d: %s vs %s (%d anchors)\n",
            tp$pair, tp$chrA, tp$chrB, tp$anchors))
st <- retain[pair == tp$pair][order(nep_chr, nep_ord)]
st[, i := .I]
st[, lost_one := xor(onA, onB)]
# null: keep exactly the same losses, randomise WHICH copy lost each one
st[, `:=`(nA = onA, nB = onB)]
fl <- st$lost_one & runif(nrow(st)) < 0.5
st[fl, `:=`(nA = onB, nB = onA)]
strip <- rbindlist(list(
  data.table(i=st$i, panel="observed",            side="A", kept=st$onA),
  data.table(i=st$i, panel="observed",            side="B", kept=st$onB),
  data.table(i=st$i, panel="same losses, shuffled between copies", side="A", kept=st$nA),
  data.table(i=st$i, panel="same losses, shuffled between copies", side="B", kept=st$nB)))
strip[, panel := factor(panel, levels=c("observed","same losses, shuffled between copies"))]
cat(sprintf("  observed: A keeps %d, B keeps %d | null: A %d, B %d\n",
            sum(st$onA), sum(st$onB), sum(st$nA), sum(st$nB)))

p2a <- ggplot(strip[kept == TRUE], aes(i, side, fill = side)) +
  geom_tile(height = 0.72) + facet_wrap(~panel, ncol = 1) +
  scale_fill_manual(values = c(A = A_COL, B = B_COL), guide = "none") +
  scale_y_discrete(limits = c("B","A")) +
  labs(title = sprintf("Gene retention along one Dionaea chromosome pair (%s / %s)", tp$chrA, tp$chrB),
       subtitle = "Each column is one ancestral gene, in Nepenthes order. A bar means that copy still has it.",
       x = "anchor genes, in order", y = NULL) + th +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        strip.text = element_text(hjust = 0, size = 10))
ggsave(file.path(OUT,"PF2a_retention_strip.pdf"), p2a, width=10, height=3.6)
ggsave(file.path(OUT,"PF2a_retention_strip.png"), p2a, width=10, height=3.6, dpi=200)

lim <- range(c(ps$more, ps$less)); lim <- c(0, lim[2]*1.05)
p2b <- ggplot(ps, aes(more, less)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(size = 3.2, colour = A_COL, alpha = 0.85) +
  annotate("text", x = lim[2]*0.55, y = lim[2]*0.72, label = "equal retention",
           colour = "grey45", size = 3.4, hjust = 0) +
  coord_fixed(xlim = lim, ylim = lim) +
  labs(title = "One subgenome consistently keeps more genes",
       subtitle = sprintf("One point per homeologous chromosome pair. %d of %d fall below the line (sign test p = %.3g).",
                          sum(ps$kept_A > ps$kept_B), nrow(ps), bt$p.value),
       x = "genes retained, more-retained copy", y = "genes retained, less-retained copy") + th
ggsave(file.path(OUT,"PF2b_fractionation_scatter.pdf"), p2b, width=6.2, height=6)
ggsave(file.path(OUT,"PF2b_fractionation_scatter.png"), p2b, width=6.2, height=6, dpi=200)
fwrite(ps, file.path(OUT,"PF2_fractionation_by_pair.csv"))
cat("  wrote PF2a, PF2b\n")

# ========================= PF3: SHUFFLE CONTROL =============================
hr("PF3. SHUFFLE CONTROL")
prop[, `:=`(g = get(gcol), lab = get(lcol), reg = get(rcol2), pos = get(mcol))]
prop <- prop[lab %in% c("A","B") & !is.na(pos)]
prop[, key := paste(genome, chr)]
runmax <- function(v) max(rle(v)$lengths)
NPERM <- 1000
sh <- prop[, {
  v <- lab[order(pos)]
  obs <- runmax(v)
  nul <- replicate(NPERM, runmax(sample(v)))
  .(n = .N, obs = obs, null_med = median(nul), null_max = max(nul),
    ratio = obs/median(nul), p = (1 + sum(nul >= obs))/(NPERM + 1))
}, by = key]
sh <- sh[n >= 100][order(-ratio)]
print(head(sh, 12), row.names=FALSE)
cat(sprintf("\n  chromosomes tested: %d | all with p <= %.4f: %s\n",
            nrow(sh), 1/(NPERM+1), all(sh$p <= 1/(NPERM+1))))
cat(sprintf("  observed longest run is %.0fx the shuffled median (median across chromosomes)\n",
            median(sh$ratio)))

TK <- "Drosera_capensis chr15_collapsed"
if (!TK %in% prop$key) { TK <- sh$key[1]; cat("  capensis chr15 absent; using", TK, "\n") }
z <- prop[key == TK][order(pos)]; z[, i := .I]
z[, shuf := sample(lab)]
bars <- rbindlist(list(
  data.table(i=z$i, panel="real labels, genes in order",    lab=z$lab),
  data.table(i=z$i, panel="same labels, shuffled positions", lab=z$shuf)))
bars[, panel := factor(panel, levels=c("real labels, genes in order","same labels, shuffled positions"))]
rr <- sh[key == TK]
p3 <- ggplot(bars, aes(i, 1, fill = lab)) +
  geom_tile() + facet_wrap(~panel, ncol = 1) +
  scale_fill_manual(values = c(A = A_COL, B = B_COL),
                    labels = c(A = "A subgenome", B = "B subgenome")) +
  labs(title = sprintf("Subgenome labels are spatially structured (%s)", sub(" ", " ", TK)),
       subtitle = sprintf("Identical A:B ratio in both rows; only the ordering differs. Longest run: %d genes observed vs %.0f shuffled.",
                          rr$obs, rr$null_med),
       x = "genes, in chromosomal order", y = NULL) + th +
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        strip.text = element_text(hjust = 0, size = 10))
ggsave(file.path(OUT,"PF3a_shuffle_strip.pdf"), p3, width=10, height=3.2)
ggsave(file.path(OUT,"PF3a_shuffle_strip.png"), p3, width=10, height=3.2, dpi=200)

sh[, sp := sub(" .*$", "", key)]
p3b <- ggplot(sh, aes(reorder(key, ratio), ratio, colour = sp)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_point(size = 2.4) + coord_flip() + scale_y_log10() +
  labs(title = "Every chromosome, observed longest run vs shuffled",
       subtitle = sprintf("%d permutations per chromosome. A ratio of 1 is what chance alone produces.", NPERM),
       x = NULL, y = "observed longest run / shuffled median") + th +
  theme(axis.text.y = element_text(size = 6))
ggsave(file.path(OUT,"PF3b_runlength_all.pdf"), p3b, width=7.5, height=9)
ggsave(file.path(OUT,"PF3b_runlength_all.png"), p3b, width=7.5, height=9, dpi=200)
fwrite(sh, file.path(OUT,"PF3_runlength.csv"))
cat("  wrote PF3a, PF3b\n")

hr("DONE")
print(list.files(OUT))
