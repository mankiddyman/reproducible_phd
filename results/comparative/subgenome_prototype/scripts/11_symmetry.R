#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); SET <- "cds"; QC <- file.path(BASE,"trees","qc")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
chk <- read_csv(file.path(QC,paste0("segment_checks_",SET,".csv")), show_col_types=FALSE)
Bl  <- read_csv(file.path(QC,paste0("panelB_blocks_",SET,".csv")), show_col_types=FALSE) %>%
       mutate(vx=as.numeric(votes_X))
G <- mean(Bl$vx)

hr("1. THE SYMMETRY TEST")
cat(sprintf("  global = %.3f | frame flip predicts X segs ~ %.3f | selection predicts ~ %.3f\n\n",
            G, G, G + (G - median(chk$other_fx[chk$call=="Y"], na.rm=TRUE))))
print(as.data.frame(chk %>% filter(!is.na(other_fx)) %>% group_by(call) %>%
  summarise(n=n(), med_other=round(median(other_fx),3),
            dev=round(median(other_fx)-G,3), .groups="drop")), row.names=FALSE)
xs <- chk$other_fx[chk$call=="X" & !is.na(chk$other_fx)]
ys <- chk$other_fx[chk$call=="Y" & !is.na(chk$other_fx)]
cat(sprintf("\n  X-seg deviation %+.3f vs Y-seg deviation %+.3f\n",
            median(xs)-G, median(ys)-G))
cat(sprintf("  |ratio| = %.2f   (near 1.0 = selection artefact; near 0 = frame flip)\n",
            abs((median(xs)-G)/(median(ys)-G))))
cat(sprintf("  Wilcoxon X segs vs global: p = %.3g\n",
            suppressWarnings(wilcox.test(xs, mu=G)$p.value)))

hr("2. is there real anchor-level shared variance?")
an <- Bl %>% group_by(anchor, region) %>%
  summarise(k=sum(vx), n=n(), fx=k/n, .groups="drop") %>% filter(n>=4)
obs_v <- var(an$fx); exp_v <- mean(G*(1-G)/an$n)
cat(sprintf("  anchors with >=4 Drosera tips: %d\n", nrow(an)))
cat(sprintf("  observed var(fx) = %.4f   binomial expectation = %.4f   ratio = %.2f\n",
            obs_v, exp_v, obs_v/exp_v))
cat("  ratio >> 1 confirms anchors differ beyond sampling — the selection\n")
cat("  artefact is then guaranteed, independent of any HE.\n\n")
cat("  per-region baselines (segments must be judged against these, not global):\n")
print(as.data.frame(an %>% group_by(region) %>%
  summarise(anchors=n(), med_fx=round(median(fx),3), mean_fx=round(mean(fx),3),
            .groups="drop") %>% arrange(mean_fx)), row.names=FALSE)

hr("3. is low fx just poor tree resolution?")
tq <- read_csv(file.path(QC,paste0("tree_qc_",SET,".csv")), show_col_types=FALSE) %>%
      select(anchor, sup_min, n_ingroup, balance)
an2 <- an %>% left_join(tq, by="anchor") %>% filter(!is.na(sup_min)) %>%
       mutate(bin=cut(fx, c(-.01,.34,.66,1.01), labels=c("Y-leaning","mixed","X-leaning")))
print(as.data.frame(an2 %>% group_by(bin) %>%
  summarise(anchors=n(), med_sup=median(sup_min), med_ntips=median(n_ingroup),
            med_balance=round(median(balance),3), .groups="drop")), row.names=FALSE)
cat(sprintf("\n  Spearman fx ~ sup_min: rho = %.3f, p = %.3g\n",
  cor(an2$fx, an2$sup_min, method="spearman"),
  suppressWarnings(cor.test(an2$fx, an2$sup_min, method="spearman")$p.value)))
cat("  Strong positive rho = Y-leaning anchors are just unresolved trees.\n")

hr("4. segment counts, region-matched baseline, SYMMETRIC filter")
rb <- an %>% group_by(region) %>% summarise(base=mean(fx), .groups="drop")
cat("  Y segment kept if other_fx is NOT depressed below its region baseline;\n")
cat("  X segment kept on the mirrored rule. Same criterion both directions.\n\n")
sym <- chk %>% filter(n_reg==1) %>%
  left_join(rb, by=c("regs"="region")) %>%
  mutate(base = coalesce(base, G),
         keep = case_when(is.na(other_fx) ~ TRUE,
                          call=="Y" ~ other_fx > base - 0.10,
                          call=="X" ~ other_fx < base + 0.10))
cat(sprintf("  single-region segments: %d   kept: %d\n\n", nrow(sym), sum(sym$keep)))
print(as.data.frame(sym %>% filter(keep) %>% count(genome, call) %>%
  pivot_wider(names_from=call, values_from=n, values_fill=0)), row.names=FALSE)
write_csv(sym, file.path(QC, paste0("segments_symmetric_", SET, ".csv")))
