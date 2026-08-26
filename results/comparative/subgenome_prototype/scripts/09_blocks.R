#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2)})
BASE <- Sys.getenv("SUBG_BASE", getwd())
SET  <- commandArgs(TRUE)[1]; if (is.na(SET)) SET <- "cds"
QC   <- file.path(BASE,"trees","qc"); NPERM <- 999; W <- 11
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))

Bl <- read_csv(file.path(QC, paste0("panelB_blocks_", SET, ".csv")), show_col_types=FALSE) %>%
  mutate(vx = as.numeric(votes_X)) %>% arrange(genome, chr_gs, ord)

## FULL windows only — no truncated edge windows
rollfull <- function(x, w) { n <- length(x); if (n < w) return(numeric(0))
  vapply(seq_len(n-w+1), function(i) mean(x[i:(i+w-1)]), numeric(1)) }
maxrun <- function(x, val) { r <- rle(x); m <- r$lengths[r$values == val]
  if (length(m)) max(m) else 0L }
stat1 <- function(v) { fr <- rollfull(v, W)
  c(run_Y = maxrun(v,0), run_X = maxrun(v,1), n_low = sum(fr <= 0.1),
    ac1 = if (sd(v) > 0) cor(v[-length(v)], v[-1]) else 0) }

hr("1. within-chromosome permutation — are Y votes clustered or scattered?")
cat(sprintf("  window %d (full only), %d perms. Null keeps each chromosome's\n", W, NPERM))
cat("  X-share fixed and destroys only the spatial ordering.\n\n")
set.seed(1)
res <- Bl %>% group_by(genome, chr_gs) %>% filter(n() >= 30) %>%
  group_modify(function(d, k) {
    v <- d$vx; obs <- stat1(v)
    nul <- vapply(seq_len(NPERM), function(i) stat1(sample(v)), numeric(4))
    tibble(n = length(v), fx = round(mean(v),3),
           run_Y = obs["run_Y"], run_Y_null = round(mean(nul["run_Y",]),1),
           p_runY = (1 + sum(nul["run_Y",] >= obs["run_Y"]))/(NPERM+1),
           n_low = obs["n_low"], p_low = (1 + sum(nul["n_low",] >= obs["n_low"]))/(NPERM+1),
           ac1 = round(obs["ac1"],3), p_ac = (1 + sum(nul["ac1",] >= obs["ac1"]))/(NPERM+1))
  }) %>% ungroup() %>% mutate(across(starts_with("p_"), ~p.adjust(.x, "BH")))

print(as.data.frame(res %>% arrange(p_runY) %>% head(20)), row.names=FALSE)
cat(sprintf("\n  chromosomes tested: %d\n", nrow(res)))
cat(sprintf("  p_runY < 0.05 (BH): %d    p_ac < 0.05: %d    p_low < 0.05: %d\n",
            sum(res$p_runY<0.05), sum(res$p_ac<0.05), sum(res$p_low<0.05)))
cat("\n  These are the chromosomes where blocks are real, not edge artefacts.\n")
write_csv(res, file.path(QC, paste0("block_permutation_", SET, ".csv")))

hr("2. per species — how many chromosomes carry contiguous Y tracts?")
print(as.data.frame(res %>% group_by(genome) %>%
  summarise(chrs=n(), sig_runY=sum(p_runY<0.05), sig_ac=sum(p_ac<0.05),
            med_fx=round(median(fx),3), .groups="drop")), row.names=FALSE)

hr("3. segment calls on the chromosomes that passed")
keep <- res %>% filter(p_runY < 0.05 | p_ac < 0.05) %>% select(genome, chr_gs)
seg <- Bl %>% semi_join(keep, by=c("genome","chr_gs")) %>%
  group_by(genome, chr_gs) %>%
  mutate(fr = c(rep(NA_real_, (W-1)%/%2), rollfull(vx, W), rep(NA_real_, W%/%2)),
         call = case_when(fr <= 0.35 ~ "Y", fr >= 0.65 ~ "X", TRUE ~ "amb")) %>%
  filter(!is.na(fr)) %>%
  mutate(seg = cumsum(call != lag(call, default = first(call)))) %>%
  group_by(genome, chr_gs, seg, call) %>%
  summarise(len = n(), start = min(ord), end = max(ord), .groups="drop") %>%
  filter(len >= 5)
cat(sprintf("  segments >= 5 anchors: %d\n\n", nrow(seg)))
print(as.data.frame(seg %>% count(genome, call) %>% pivot_wider(names_from=call,
      values_from=n, values_fill=0)), row.names=FALSE)
cat("\n  Both X and Y segments in a species = both lineages present, in blocks.\n")
cat("  This is the block-level version of the claim, replacing chromosome labels.\n")
write_csv(seg, file.path(QC, paste0("segments_", SET, ".csv")))
