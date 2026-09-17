#!/usr/bin/env Rscript
# ============================================================================
# DR06 — WHAT DO THE 1,232 GENE TREES ACTUALLY SAY?
#
# Support values saturate (UFBoot 100/100 everywhere, local PP 1.0 everywhere)
# because there are 2 Mb of sites and 1,232 gene trees. Gene-tree COUNTS do
# not saturate. This script reads them directly.
#
# THE TEST THAT MATTERS — MINORITY SYMMETRY
#   Under pure incomplete lineage sorting the coalescent is SYMMETRIC: when
#   lineages fail to coalesce, the two wrong resolutions are equally likely.
#   So the two MINORITY counts must be EQUAL.
#     symmetric   -> ILS. Biology. Not fixable, but modellable (ASTRAL).
#     ASYMMETRIC  -> something non-coalescent: PARALOGY (mixing A1 and A2 into
#                    one tip), introgression, or systematic error. FIXABLE.
#   With ~600 minority trees a 5-point asymmetry is detectable at 80% power.
#   This is a free, hard test and it has never been run on this dataset.
#
# ALSO: coalescent branch length predicts concordance exactly.
#     P(gene tree matches) = 1 - (2/3)exp(-t),  t in coalescent units
#   ASTRAL gave regia's branch at 0.198-0.250 CU -> 45-48% expected.
#   If OBSERVED concordance is far from that, the MSC does not fit and the
#   conflict is not ILS.
#
# IN   DR/tree/genetrees/*.treefile
# OUT  DR/out/DR06_quartet_counts.csv, DR06_regia.csv
# FIG  DR/fig/DR06_1_quartets.pdf, DR06_2_regia.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(ape); library(dplyr); library(readr); library(tidyr); library(ggplot2)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

hr("0. LOAD")
fs <- list.files("DR/tree/genetrees", pattern="[.]treefile$", full.names=TRUE)
cat(sprintf("  gene tree files: %d\n", length(fs)))
TR <- lapply(fs, function(f) tryCatch(read.tree(f), error=function(e) NULL))
names(TR) <- sub("[.]treefile$", "", basename(fs))
TR <- TR[!vapply(TR, is.null, logical(1))]
cat(sprintf("  parsed: %d\n", length(TR)))
cat("  tips per gene tree:\n")
print(summary(vapply(TR, function(t) length(t$tip.label), integer(1))))

# quartet resolution for four named tips in one unrooted gene tree
quart <- function(t, a, b, c, d) {
  tips <- c(a,b,c,d)
  if (!all(tips %in% t$tip.label)) return(NA_character_)
  q <- drop.tip(t, setdiff(t$tip.label, tips))
  if (length(q$tip.label) != 4) return(NA_character_)
  # for an unrooted 4-tip tree, the internal split defines the pairing
  sp <- prop.part(unroot(q))
  s  <- q$tip.label[sp[[length(sp)]]]
  if (length(s) != 2) return("unresolved")
  paste(sort(s), collapse="|")
}

count_quartet <- function(a, b, c, d, label) {
  v <- vapply(TR, function(t) quart(t, a, b, c, d), character(1))
  v <- v[!is.na(v) & v != "unresolved"]
  if (!length(v)) return(NULL)
  tb <- sort(table(v), decreasing=TRUE)
  n  <- sum(tb)
  r  <- as.integer(tb); r <- c(r, rep(0L, max(0, 3-length(r))))[1:3]
  # symmetry test on the two minorities: hard prediction of the coalescent
  bt <- if (r[2]+r[3] >= 10) binom.test(r[2], r[2]+r[3], 0.5)$p.value else NA_real_
  tibble(test=label, n=n,
         top=names(tb)[1], top_n=r[1], top_f=r[1]/n,
         mid=ifelse(length(tb)>1, names(tb)[2], NA), mid_n=r[2],
         low=ifelse(length(tb)>2, names(tb)[3], NA), low_n=r[3],
         minority_p=bt,
         ILS_consistent=ifelse(is.na(bt), NA, bt > 0.05))
}

hr("1. REGIA — the contested node, both subgenomes, all core species")
CORE <- c("binata","paradoxa","scorpioides","capensis")
RG <- bind_rows(lapply(c("A","B"), function(s)
  bind_rows(lapply(CORE, function(cs)
    count_quartet(paste0("regia_",s), paste0("Dionaea_",s),
                  paste0(cs,"_",s), "Nepenthes",
                  sprintf("regia_%s vs %s_%s", s, cs, s))))))
RG <- RG %>% mutate(
  regia_Dionaea = grepl("Dionaea", top) & grepl("regia", top),
  across(c(top_f), ~round(.,3)), minority_p=signif(minority_p,3))
print(as.data.frame(RG %>% select(test, n, top, top_f, mid_n, low_n,
                                  minority_p, ILS_consistent)), row.names=FALSE)
cat("\n  top = the pairing recovered most often. top_f = its frequency.\n")
cat("  minority_p tests whether the two LOSING resolutions are equally common.\n")
cat("  p > 0.05 -> symmetric -> consistent with pure ILS.\n")
cat("  p < 0.05 -> ASYMMETRIC -> paralogy, introgression or systematic error.\n")
nsym <- sum(!RG$ILS_consistent, na.rm=TRUE)
cat(sprintf("\n  asymmetric quartets: %d of %d\n", nsym, sum(!is.na(RG$ILS_consistent))))

hr("2. POOLED, per subgenome")
pool <- RG %>% mutate(sub=substr(sub("^regia_", "", test),1,1)) %>%
  group_by(sub) %>%
  summarise(quartets=n(), mean_top_f=round(mean(top_f),3),
            regia_Dionaea_wins=sum(regia_Dionaea),
            total_minority=sum(mid_n+low_n),
            .groups="drop")
print(as.data.frame(pool), row.names=FALSE)

hr("3. DOES THE OBSERVED CONCORDANCE MATCH THE COALESCENT PREDICTION?")
cat("  P(match) = 1 - (2/3)exp(-t). ASTRAL put regia's branch at 0.198 (A)\n")
cat("  and 0.250 (B) coalescent units.\n\n")
for (x in list(c("A",0.198), c("B",0.250))) {
  t <- as.numeric(x[2]); pred <- 1 - (2/3)*exp(-t)
  obs <- mean(RG$top_f[substr(sub("^regia_","",RG$test),1,1)==x[1]])
  cat(sprintf("  subgenome %s: %.3f CU -> predicted %.3f | observed %.3f | diff %+.3f\n",
              x[1], t, pred, obs, obs-pred))
}
cat("\n  close agreement -> the MSC fits and the conflict IS ILS.\n")
cat("  observed much LOWER -> more conflict than ILS explains.\n")

hr("4. CONTROL — nodes we are confident about")
cat("  If the method works, uncontested nodes should show HIGH top_f and\n")
cat("  symmetric minorities. If they do not, the gene trees are just noisy.\n\n")
CTL <- bind_rows(
  count_quartet("paradoxa_A","scorpioides_A","binata_A","Nepenthes",
                "(paradoxa,scorpioides) vs binata  [A]"),
  count_quartet("paradoxa_B","scorpioides_B","binata_B","Nepenthes",
                "(paradoxa,scorpioides) vs binata  [B]"),
  count_quartet("binata_A","paradoxa_A","capensis_A","Nepenthes",
                "core vs capensis  [A]"),
  count_quartet("Dionaea_A","Dionaea_B","binata_A","Nepenthes",
                "Dionaea A|B split  [the allopolyploidy]"))
print(as.data.frame(CTL %>% mutate(top_f=round(top_f,3),
  minority_p=signif(minority_p,3)) %>%
  select(test, n, top_f, mid_n, low_n, minority_p, ILS_consistent)),
  row.names=FALSE)

hr("5. ALL THREE RESOLUTIONS, side by side")
LONG <- RG %>%
  transmute(test, `1st`=top_n, `2nd`=mid_n, `3rd`=low_n) %>%
  pivot_longer(-test, names_to="rank", values_to="n")
p1 <- ggplot(LONG, aes(reorder(test, n), n, fill=rank)) +
  geom_col(position="fill", width=0.7) +
  geom_hline(yintercept=c(1/3,2/3), linetype="dotted", colour="grey30") +
  coord_flip() +
  scale_fill_manual(values=c("1st"="#1D9E75","2nd"="grey65","3rd"="grey80")) +
  labs(title="DR06 - gene-tree quartet counts at regia's node",
       subtitle="dotted = 1/3 and 2/3 | equal 2nd and 3rd bars = consistent with ILS",
       x=NULL, y="fraction of gene trees", fill="resolution") +
  theme_minimal(9) + theme(legend.position="top")
suppressWarnings(ggsave("DR/fig/DR06_1_quartets.pdf", p1, width=9, height=6))
write_csv(bind_rows(RG %>% mutate(set="regia"), CTL %>% mutate(set="control")),
          "DR/out/DR06_quartet_counts.csv")

hr("VERDICT")
mt <- mean(RG$top_f)
cat(sprintf("  mean top resolution frequency at regia's node: %.3f\n", mt))
cat(sprintf("  regia+Dionaea is the top resolution in %d of %d quartets\n",
            sum(RG$regia_Dionaea), nrow(RG)))
cat(sprintf("  minorities symmetric in %d of %d\n",
            sum(RG$ILS_consistent, na.rm=TRUE), sum(!is.na(RG$ILS_consistent))))
cat("\n  READ:\n")
cat("   top_f near 0.33  -> no signal; the node may be a genuine polytomy\n")
cat("   top_f 0.45-0.55 with SYMMETRIC minorities -> real branch, heavy ILS.\n")
cat("     Fixable only by modelling (ASTRAL), which is what we did.\n")
cat("   top_f 0.45-0.55 with ASYMMETRIC minorities -> NOT pure ILS. Paralogy\n")
cat("     from mixing A1 and A2 into one tip is the leading candidate, and\n")
cat("     that IS fixable -- step 3 tests it directly.\n")
