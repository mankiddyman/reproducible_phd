#!/usr/bin/env Rscript
# Stage 03 — the gate.  Rscript 03_treeqc.R {cds|aa|p12}
# Decisive number: UFBoot on the DEEPEST INGROUP SPLIT.
# ape::root() shuffles node labels, so support is read off the tree AS STORED,
# keyed by bipartition, then looked up again after rooting.
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({ library(ape); library(dplyr); library(readr); library(ggplot2) })

BASE <- Sys.getenv("SUBG_BASE", getwd())
SET  <- commandArgs(TRUE)[1]; if (is.na(SET)) SET <- "cds"
TRE  <- file.path(BASE, "trees", paste0("tre_", SET))
OUTD <- file.path(BASE, "trees", "qc"); dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
GENOMES <- c("Dionaea_muscipula","Drosera_binata","Drosera_capensis",
             "Drosera_paradoxa","Drosera_regia","Drosera_scorpioides","Nepenthes_gracilis")
NEP <- "Nepenthes_gracilis"

genome_of <- function(tips) {
  hit <- rep(NA_character_, length(tips))
  for (g in GENOMES) hit[is.na(hit) & startsWith(tips, g)] <- g
  hit
}
split_key <- function(tipset, all_tips) {
  a <- sort(tipset); b <- sort(setdiff(all_tips, tipset))
  ka <- paste(a, collapse="|"); kb <- paste(b, collapse="|"); if (ka < kb) ka else kb
}
support_map <- function(tr) {
  lab <- suppressWarnings(as.numeric(tr$node.label)); desc <- prop.part(tr); all_t <- tr$tip.label
  m <- new.env(hash=TRUE, parent=emptyenv())
  for (i in seq_along(desc)) if (!is.na(lab[i]))
    assign(split_key(all_t[desc[[i]]], all_t), lab[i], envir=m)
  m
}
lookup <- function(m, tipset, all_tips) {
  k <- split_key(tipset, all_tips)
  if (exists(k, envir=m, inherits=FALSE)) get(k, envir=m) else NA_real_
}

files <- list.files(TRE, pattern="\\.treefile$", full.names=TRUE)
cat(sprintf("set = %s   treefiles = %d\n", SET, length(files)))
if (!length(files)) stop("no treefiles")
rows <- vector("list", length(files)); tiprows <- vector("list", length(files))

for (i in seq_along(files)) {
  anchor <- sub("\\.treefile$", "", basename(files[i]))
  tr <- tryCatch(read.tree(files[i]), error=function(e) NULL)
  if (is.null(tr) || Ntip(tr) < 4) next
  smap <- support_map(tr); all_t <- tr$tip.label; gen <- genome_of(all_t)
  nept <- all_t[gen == NEP]; if (length(nept) != 1) next
  rt <- tryCatch(root(tr, outgroup=nept, resolve.root=TRUE), error=function(e) NULL)
  if (is.null(rt)) next
  rn <- Ntip(rt) + 1L; kids <- rt$edge[rt$edge[,1] == rn, 2]
  ing <- kids[kids != match(nept, rt$tip.label)]
  if (length(ing) != 1 || ing <= Ntip(rt)) next
  sub_ <- extract.clade(rt, ing); kid2 <- rt$edge[rt$edge[,1] == ing, 2]
  if (length(kid2) != 2) next
  side_tips <- lapply(kid2, function(k)
    if (k <= Ntip(rt)) rt$tip.label[k] else extract.clade(rt, k)$tip.label)
  L <- side_tips[[1]]; R <- side_tips[[2]]
  supL <- lookup(smap, L, all_t); supR <- lookup(smap, R, all_t)
  nep_bl <- rt$edge.length[rt$edge[,2] == match(nept, rt$tip.label)]
  depth <- max(node.depth.edgelength(sub_))
  r2t <- setNames(node.depth.edgelength(rt)[seq_len(Ntip(rt))], rt$tip.label)
  ing_g <- gen[match(c(L,R), all_t)]; n_ing <- length(L) + length(R)
  rows[[i]] <- tibble(anchor, n_ingroup=n_ing, a=length(L), b=length(R),
    balance=min(length(L),length(R))/n_ing, sup_L=supL, sup_R=supR,
    sup_min=suppressWarnings(min(supL,supR,na.rm=TRUE)),
    nep_bl=if (length(nep_bl)) nep_bl else NA_real_, ingroup_depth=depth,
    nep_ratio=if (length(nep_bl) && depth>0) nep_bl/depth else NA_real_,
    n_species=length(unique(ing_g)))
  tiprows[[i]] <- tibble(anchor, tip=c(L,R), genome=ing_g,
    side=rep(c("L","R"), c(length(L),length(R))), r2t=unname(r2t[c(L,R)]),
    n_ingroup=n_ing, a=length(L), b=length(R))
}
trees <- bind_rows(rows); tips <- bind_rows(tiprows)
trees$sup_min[is.infinite(trees$sup_min)] <- NA_real_

sp <- tips %>% group_by(anchor, genome) %>%
  summarise(k=n(), in_L=sum(side=="L"), .groups="drop") %>%
  left_join(trees %>% select(anchor,n_ingroup,a,b), by="anchor") %>%
  mutate(spans = ifelse(k>=2, in_L>0 & in_L<k, NA),
         p_null = ifelse(k>=2 & k<=n_ingroup,
                         1-(choose(a,k)+choose(b,k))/choose(n_ingroup,k), NA_real_))

write_csv(trees, file.path(OUTD, paste0("tree_qc_", SET, ".csv")))
write_csv(tips,  file.path(OUTD, paste0("tip_sides_", SET, ".csv")))
write_csv(sp,    file.path(OUTD, paste0("species_spanning_", SET, ".csv")))

hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
hr("1. THE GATE — UFBoot on the deepest ingroup split")
q <- quantile(trees$sup_min, c(.1,.25,.5,.75,.9), na.rm=TRUE)
cat(sprintf("  trees analysed : %d of %d\n", nrow(trees), length(files)))
cat(sprintf("  support NA     : %d\n", sum(is.na(trees$sup_min))))
cat("  quantiles      : "); print(round(q,1))
cat(sprintf("  >= 95 : %.3f    >= 80 : %.3f    < 70 : %.3f\n",
  mean(trees$sup_min>=95,na.rm=TRUE), mean(trees$sup_min>=80,na.rm=TRUE),
  mean(trees$sup_min<70,na.rm=TRUE)))
med <- q[["50%"]]
cat(sprintf("\n  VERDICT: median %.0f -> %s\n", med,
  if (med>=80) "PROCEED. Deep split is resolved."
  else if (med>=70) "MARGINAL. Filter to sup_min >= 80 downstream and report n."
  else "STOP. Deep node unresolved. Check AA before abandoning."))

hr("2. rooting sanity")
cat(sprintf("  Nepenthes branch / ingroup depth: median %.2f  (>2 = long-branch outgroup)\n",
  median(trees$nep_ratio, na.rm=TRUE)))
cat(sprintf("  trees with ratio > 2: %d (%.1f%%)\n", sum(trees$nep_ratio>2,na.rm=TRUE),
  100*mean(trees$nep_ratio>2,na.rm=TRUE)))
cat("  NOTE: one Nepenthes tip means outgroup rooting cannot fail. Plausibility\n")
cat("        check only; midpoint comparison is a separate step.\n")

hr("3. tree shape census — this is what sets the null")
cat("  ingroup tips (n):\n"); print(summary(trees$n_ingroup))
cat("\n  root split balance min(a,b)/n:\n"); print(round(summary(trees$balance),3))
cat(sprintf("\n  very lopsided (balance < 0.2): %.1f%%\n", 100*mean(trees$balance<0.2)))
cat("  species present per tree:\n"); print(table(trees$n_species))

hr("4. copy number and spanning — observed vs EXACT null")
tab <- sp %>% filter(genome != NEP) %>% group_by(genome) %>%
  summarise(trees_present=n(), mean_k=round(mean(k),2), frac_k1=round(mean(k==1),3),
    n_testable=sum(k>=2), obs_span=round(mean(spans,na.rm=TRUE),3),
    exp_span=round(mean(p_null,na.rm=TRUE),3),
    z=round((sum(spans,na.rm=TRUE)-sum(p_null,na.rm=TRUE))/
            sqrt(sum(p_null*(1-p_null),na.rm=TRUE)),2), .groups="drop") %>%
  arrange(desc(obs_span))
print(as.data.frame(tab), row.names=FALSE)
cat("\n  exp_span is NOT 0.5 — it is set per tree by n, the a/b split, and k.\n")
cat("  A 4-copy species spans ~93% of the time by chance alone. Read z.\n")

hr("5. rate control — is spanning just long branches?")
rate <- tips %>% group_by(anchor, genome) %>%
  summarise(k=n(), rate_gap=if (n()>=2) max(r2t)-min(r2t) else NA_real_,
            mean_r2t=mean(r2t), .groups="drop") %>%
  inner_join(sp %>% select(anchor,genome,spans), by=c("anchor","genome")) %>%
  filter(!is.na(spans), !is.na(rate_gap))
if (nrow(rate) > 30) {
  w <- wilcox.test(rate_gap ~ spans, data=rate)
  g <- glm(spans ~ scale(rate_gap) + scale(mean_r2t), data=rate, family=binomial)
  cat(sprintf("  Wilcoxon rate_gap ~ spans : p = %.3g\n", w$p.value))
  cat("  logistic spans ~ rate_gap + mean_r2t:\n")
  print(round(summary(g)$coefficients, 4))
  cat("\n  Label permutation cannot see branch lengths, so this is the only\n")
  cat("  check on LBA. Significant rate_gap = spanning is measuring rate.\n")
} else cat("  too few testable rows\n")

p <- ggplot(trees, aes(sup_min)) +
  geom_histogram(binwidth=5, fill="#4C6EF5", colour="white") +
  geom_vline(xintercept=c(70,80,95), linetype="dashed", colour="firebrick") +
  labs(title=paste0("Deepest ingroup split — UFBoot (", SET, ")"),
       subtitle="dashed: 70 / 80 / 95. Median below 80 = gate does not pass.",
       x="UFBoot", y="trees") + theme_minimal(11)
suppressWarnings(ggsave(file.path(OUTD, paste0("QC_support_", SET, ".pdf")), p, width=7, height=4.5))
cat(sprintf("\nWROTE: %s/{tree_qc,tip_sides,species_spanning}_%s.csv  QC_support_%s.pdf\n",
  OUTD, SET, SET))
