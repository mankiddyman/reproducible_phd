#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd())
SET  <- commandArgs(TRUE)[1]; if (is.na(SET)) SET <- "cds"
QC   <- file.path(BASE,"trees","qc"); W <- 11
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
rollfull <- function(x,w){n<-length(x); if(n<w) return(rep(NA_real_,n))
  c(rep(NA_real_,(w-1)%/%2), vapply(seq_len(n-w+1),function(i)mean(x[i:(i+w-1)]),numeric(1)), rep(NA_real_,w%/%2))}

Bl  <- read_csv(file.path(QC,paste0("panelB_blocks_",SET,".csv")), show_col_types=FALSE) %>%
       mutate(vx=as.numeric(votes_X)) %>% arrange(genome, chr_gs, ord)
res <- read_csv(file.path(QC,paste0("block_permutation_",SET,".csv")), show_col_types=FALSE)
keep <- res %>% filter(p_runY < 0.05 | p_ac < 0.05) %>% select(genome, chr_gs)
GLOBAL <- mean(Bl$vx)
cat(sprintf("global X-share = %.3f   chromosomes passing = %d\n", GLOBAL, nrow(keep)))

hr("1. CONFOUND A — do call changes just track ancestral-region changes?")
cat("  scorpioides median 6 regions/chr, paradoxa 3.5 (v2 8.3). A step at a\n")
cat("  region boundary is a change of reference region, not of A/B ancestry.\n\n")
tr <- Bl %>% semi_join(keep, by=c("genome","chr_gs")) %>%
  group_by(genome, chr_gs) %>%
  mutate(fr = rollfull(vx, W),
         call = case_when(fr <= .35 ~ "Y", fr >= .65 ~ "X", TRUE ~ "amb")) %>%
  filter(!is.na(fr)) %>%
  mutate(call_chg = call != lag(call), reg_chg = region != lag(region)) %>%
  filter(!is.na(call_chg)) %>% ungroup()
tab <- table(call_change = tr$call_chg, region_change = tr$reg_chg)
print(tab)
ft <- fisher.test(tab)
cat(sprintf("\n  Fisher OR = %.2f, p = %.3g\n", ft$estimate, ft$p.value))
cat("  OR >> 1 = the segmentation is largely region structure. OR ~ 1 = independent.\n\n")
print(as.data.frame(tr %>% group_by(genome) %>%
  summarise(anchors=n(), call_chgs=sum(call_chg), reg_chgs=sum(reg_chg),
            both=sum(call_chg & reg_chg),
            frac_call_at_reg=round(sum(call_chg & reg_chg)/max(1,sum(call_chg)),3),
            .groups="drop")), row.names=FALSE)

hr("2. CONFOUND B — are Y segments Dionaea frame flips?")
cat("  At a segment's own anchors, what do the OTHER species do?\n")
cat(sprintf("  Others near %.2f = species-specific block. Others also low = frame flip.\n\n", GLOBAL))
segs <- tr %>% group_by(genome, chr_gs) %>%
  mutate(seg = cumsum(call != lag(call, default=first(call)))) %>%
  group_by(genome, chr_gs, seg, call) %>%
  filter(n() >= 5, call != "amb") %>%
  summarise(len=n(), n_reg=n_distinct(region), regs=paste(unique(sub("_dom","",region)),collapse=","),
            anchors=list(anchor), .groups="drop")
chk <- bind_rows(lapply(seq_len(nrow(segs)), function(i) {
  a <- segs$anchors[[i]]
  o <- Bl %>% filter(anchor %in% a, genome != segs$genome[i]) %>%
       group_by(genome) %>% summarise(fx=mean(vx), n=n(), .groups="drop") %>% filter(n>=3)
  segs[i,] %>% select(genome, chr_gs, call, len, n_reg, regs) %>%
    mutate(self_n=length(a), n_other_sp=nrow(o),
           other_fx=if(nrow(o)) round(mean(o$fx),3) else NA_real_,
           n_other_low=if(nrow(o)) sum(o$fx<0.5) else NA_integer_)
}))
cat("  --- Y segments ---\n")
print(as.data.frame(chk %>% filter(call=="Y") %>% arrange(other_fx)), row.names=FALSE)
yy <- chk %>% filter(call=="Y", !is.na(other_fx))
cat(sprintf("\n  Y segments: %d   median other_fx = %.3f   with >=3 other species also Y-leaning: %d\n",
            nrow(yy), median(yy$other_fx), sum(yy$n_other_low>=3, na.rm=TRUE)))
cat(sprintf("  binomial vs global %.3f: p = %.3g\n", GLOBAL,
    if(nrow(yy)) wilcox.test(yy$other_fx, mu=GLOBAL)$p.value else NA))

hr("3. surviving segments — single-region, not frame-flipped")
surv <- chk %>% filter(n_reg == 1, is.na(other_fx) | other_fx > 0.5 | call=="X")
cat(sprintf("  segments retained: %d of %d\n\n", nrow(surv), nrow(chk)))
print(as.data.frame(surv %>% count(genome, call) %>%
  pivot_wider(names_from=call, values_from=n, values_fill=0)), row.names=FALSE)
cat("\n  Both X and Y here = both lineages present in blocks, region-confound\n")
cat("  and frame-flip both excluded. This is the defensible version.\n")
write_csv(chk, file.path(QC, paste0("segment_checks_", SET, ".csv")))
