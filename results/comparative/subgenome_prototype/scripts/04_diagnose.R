#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr)})
BASE <- Sys.getenv("SUBG_BASE", getwd())
SET <- commandArgs(TRUE)[1]; if (is.na(SET)) SET <- "cds"
QC <- file.path(BASE,"trees","qc")
trees <- read_csv(file.path(QC,paste0("tree_qc_",SET,".csv")), show_col_types=FALSE)
tips  <- read_csv(file.path(QC,paste0("tip_sides_",SET,".csv")), show_col_types=FALSE)
sp    <- read_csv(file.path(QC,paste0("species_spanning_",SET,".csv")), show_col_types=FALSE)
NEP <- "Nepenthes_gracilis"
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))

zfun <- function(d) {
  m <- !is.na(d$spans) & !is.na(d$p_null); n <- sum(m)
  if (n < 10) return(tibble(n=n, obs=NA_real_, exp=NA_real_, z=NA_real_))
  o <- sum(d$spans[m]); e <- sum(d$p_null[m]); v <- sum(d$p_null[m]*(1-d$p_null[m]))
  tibble(n=n, obs=round(o/n,3), exp=round(e/n,3), z=round((o-e)/sqrt(v),2))
}

hr("1. anatomy of the deepest split")
trees <- trees %>% mutate(minside = pmin(a,b))
cat("  size of the smaller side:\n"); print(table(trees$minside))
cat(sprintf("\n  1-vs-rest splits : %.1f%%\n", 100*mean(trees$minside==1)))
cat(sprintf("  smaller side <= 2: %.1f%%\n", 100*mean(trees$minside<=2)))

hr("2. is the lone basal tip just the longest branch?")
lone <- tips %>% group_by(anchor) %>%
  mutate(nL=sum(side=="L"), nR=sum(side=="R"),
         is_lone=(side=="L" & nL==1)|(side=="R" & nR==1),
         rank_r2t=rank(-r2t), n_ing=n()) %>% ungroup() %>% filter(is_lone)
cat(sprintf("  lone tip is the LONGEST branch in its tree: %.1f%% of cases\n",
            100*mean(lone$rank_r2t==1)))
cat(sprintf("  chance expectation (1/n): %.1f%%\n", 100*mean(1/lone$n_ing)))
cat("  genome owning the lone basal tip:\n")
print(round(prop.table(table(lone$genome)),3))

hr("3. does lopsidedness track low support?")
print(as.data.frame(trees %>%
  mutate(bin=cut(balance, c(0,.15,.25,.35,.5), include.lowest=TRUE)) %>%
  group_by(bin) %>% summarise(n=n(), med_sup=median(sup_min,na.rm=TRUE),
                              med_maxr2t=NA, .groups="drop")), row.names=FALSE)

rt <- tips %>% group_by(anchor, genome) %>%
  summarise(rate_gap=if (n()>=2) max(r2t)-min(r2t) else NA_real_,
            mean_r2t=mean(r2t), .groups="drop")
d <- sp %>% inner_join(rt, by=c("anchor","genome")) %>%
  filter(genome != NEP, !is.na(spans), !is.na(rate_gap))

hr("4. spanning stratified by within-species rate gap (tertiles)")
cat("  If the signal is ancestry it survives in tertile 1 (smallest gap).\n")
cat("  If it is LBA, z collapses toward 0 there.\n\n")
print(as.data.frame(d %>% group_by(genome) %>% mutate(tert=ntile(rate_gap,3)) %>%
  group_by(genome, tert) %>% group_modify(~zfun(.x)) %>% ungroup()), row.names=FALSE)

hr("5. clean subset: sup_min>=80, balance>=0.3, all 6 species")
keep <- trees %>% filter(sup_min>=80, balance>=0.3, n_species==6) %>% pull(anchor)
cat(sprintf("  trees kept: %d of %d\n\n", length(keep), nrow(trees)))
dd <- d %>% filter(anchor %in% keep)
print(as.data.frame(dd %>% group_by(genome) %>% group_modify(~zfun(.x)) %>% ungroup()),
      row.names=FALSE)
if (nrow(dd) > 40) {
  g <- glm(spans ~ scale(rate_gap) + scale(mean_r2t), data=dd, family=binomial)
  cat("\n  rate control ON THE CLEAN SUBSET:\n"); print(round(summary(g)$coefficients,4))
  cat("\n  rate_gap ns here = filtering fixed it. Still significant = it did not.\n")
}
