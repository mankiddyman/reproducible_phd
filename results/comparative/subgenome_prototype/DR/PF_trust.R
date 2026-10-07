#!/usr/bin/env Rscript
# ============================================================================
# PF_trust -- is the per-gene Delta interleaving what the noise model predicts,
#             and do RAW votes (not propagated labels) cluster along the chromosome?
#
# 1  minority-vote fraction inside each called block vs the predicted 0.31
#    (A at -0.16, B at +0.16, per-gene sd 0.32 -> P(wrong side) = P(Z>0.5))
# 2  run-length test on the RAW sign of Delta. The earlier test used propagated
#    labels, where every gene in a block shares the block's call, so long runs
#    were guaranteed. This one uses each gene's own vote.
# 3  how many genes must be averaged before the call is stable
# ============================================================================
suppressMessages(library(data.table))
setwd(Sys.getenv("SUBG_BASE", getwd())); options(width=205)
set.seed(3)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))
NPERM <- 499

hr("0. JOIN RAW DELTA TO THE BLOCK CALLS")
G <- fread("DR/out/DR02_gene_delta.csv")[is.finite(delta)]
G[, mb := mid/1e6]
P <- fread("DR/out/DR02_propagated_blocks.csv")[label %in% c("A","B")]
key <- if ("gene" %in% names(P)) "gene" else "tip"
M <- merge(G[, .(genome, chr, gene, mb, delta)],
           P[, .(genome, chr, gene = get(key), region, label)],
           by=c("genome","chr","gene"))
cat(sprintf("  genes with BOTH a raw Delta and a block call: %s\n",
            format(nrow(M), big.mark=",")))
M[, vote := fifelse(delta < 0, "A", "B")]
M[, agree := vote == label]
cat(sprintf("  overall agreement between a gene's own vote and its block: %.3f\n",
            mean(M$agree)))

hr("1. MINORITY VOTES INSIDE BLOCKS")
pred <- pnorm(0.16/0.32, lower.tail=FALSE)
cat(sprintf("  predicted from the DR02 noise model: %.3f of genes vote against\n", pred))
BL <- M[, .(n=.N, minority=mean(!agree), med=median(delta)), by=.(genome, chr, region, label)]
BL <- BL[n >= 10]
cat(sprintf("  blocks: %d | minority fraction: median %.3f | IQR %.3f-%.3f\n",
    nrow(BL), median(BL$minority), quantile(BL$minority,.25), quantile(BL$minority,.75)))
cat(sprintf("  blocks within 0.05 of the prediction: %d (%.0f%%)\n",
    sum(abs(BL$minority - pred) < 0.05), 100*mean(abs(BL$minority - pred) < 0.05)))
cat("\n  by species:\n")
print(M[, .(genes=.N, minority=round(mean(!agree),3)), by=genome][order(minority)],
      row.names=FALSE)
cat("\n  worst 6 blocks (most internal disagreement):\n")
print(head(BL[order(-minority), .(genome, chr, region, label, n,
                                  minority=round(minority,3), med=round(med,3))], 6),
      row.names=FALSE)

hr("2. RUN-LENGTH ON RAW VOTES -- the non-circular test")
cat("  earlier panel used propagated labels: every gene in a block shares the\n")
cat("  block call, so runs were long BY CONSTRUCTION. This uses raw sign only.\n\n")
M[, ck := paste(genome, chr)]
mx <- function(v) max(rle(v)$lengths)
R <- rbindlist(lapply(unique(M$ck), function(k){
  v <- M[ck==k][order(mb), vote]
  if (length(v) < 50 || uniqueN(v) < 2) return(NULL)
  o <- mx(v); nul <- replicate(NPERM, mx(sample(v)))
  data.table(chrom=k, n=length(v), obs=o, null_med=median(nul),
             null_p95=as.numeric(quantile(nul,.95)),
             p=(1+sum(nul>=o))/(NPERM+1)) }))
R[, ratio := obs/pmax(null_med,1)]
setorder(R, -ratio)
cat(sprintf("  chromosomes tested: %d\n", nrow(R)))
cat(sprintf("  real longest run of identical raw votes: median %d | shuffled %d\n",
    as.integer(median(R$obs)), as.integer(median(R$null_med))))
cat(sprintf("  ratio real:shuffled -- median %.2fx | range %.2fx-%.2fx\n",
    median(R$ratio), min(R$ratio), max(R$ratio)))
cat(sprintf("  chromosomes with p <= 0.05: %d of %d\n", sum(R$p <= 0.05), nrow(R)))
print(head(R, 6), row.names=FALSE)
cat("\n  NOTE: raw votes are ~31%% wrong at random, so even a perfect block gives\n")
cat("  short runs. A modest ratio here is expected; the signal lives in the MEAN,\n")
cat("  which section 3 measures directly.\n")

hr("3. HOW MANY GENES BEFORE THE CALL IS STABLE")
for (w in c(1, 5, 10, 25, 50, 100, 200)) {
  acc <- M[, { v <- .SD[order(mb)]
    if (.N < w) NULL else {
      k <- floor(.N/w)
      idx <- rep(seq_len(k), each=w)[seq_len(k*w)]
      d <- data.table(dl = v$delta[seq_len(k*w)], lb = v$label[seq_len(k*w)], g = idx)
      s <- d[, .(m = median(dl), lb = names(sort(table(lb), decreasing=TRUE))[1]), by=g]
      s[, .(ok = (m < 0) == (lb == "A"))] } }, by=ck]
  cat(sprintf("  windows of %3d genes: call matches the block %.1f%% of the time\n",
      w, 100*mean(acc$ok))) }
cat("\n  a single gene is near-useless; a block of hundreds is not.\n")
