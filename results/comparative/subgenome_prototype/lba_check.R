#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Does the Dionaea split track ANCESTRY or just BRANCH LENGTH?
# If separation happens mainly when the two copies differ most in root-to-tip
# distance, it is long-branch attraction, not subgenome signal.
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ape) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")

ref <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE) %>%
  filter(anchor == "Nepenthes_gracilis")
xy <- bind_rows(transmute(ref, chr = winner, side = "X"),
                transmute(ref, chr = ifelse(winner == chrA, chrB, chrA), side = "Y")) %>%
  mutate(tag = paste0("Dio_", sub("^chr([0-9]+)_", "c\\1_", chr)))

d <- bind_rows(lapply(list.files("win25/tre", "\\.treefile$", full.names = TRUE), function(f) {
  tr <- tryCatch(read.tree(f), error = function(e) NULL)
  if (is.null(tr) || !("NEP" %in% tr$tip.label)) return(NULL)
  tr <- root(tr, "NEP", resolve.root = TRUE)
  dx <- intersect(tr$tip.label, xy$tag[xy$side == "X"])
  dy <- intersect(tr$tip.label, xy$tag[xy$side == "Y"])
  if (length(dx) != 1 || length(dy) != 1) return(NULL)
  ing <- drop.tip(tr, "NEP")
  k <- ing$edge[ing$edge[, 1] == Ntip(ing) + 1, 2]
  if (length(k) != 2) return(NULL)
  h1 <- if (k[1] <= Ntip(ing)) ing$tip.label[k[1]] else extract.clade(ing, k[1])$tip.label
  rt <- setNames(diag(vcv(tr)), tr$tip.label)          # root-to-tip distance
  ot <- setdiff(tr$tip.label, c("NEP", dx, dy))
  tibble(win = sub("\\.treefile$", "", basename(f)),
         sep = (dx %in% h1) != (dy %in% h1),
         dlen = abs(rt[dx] - rt[dy]),
         rel  = abs(rt[dx] - rt[dy]) / mean(c(rt[dx], rt[dy])),
         dio_long = as.numeric(mean(c(rt[dx], rt[dy])) > median(rt[ot])))
}))

cat(sprintf("windows: %d | separated: %d\n", nrow(d), sum(d$sep)))
cat("\n=== root-to-tip difference between the two Dionaea copies ===\n")
print(as.data.frame(d %>% group_by(sep) %>%
  summarise(n = n(), med_abs = median(dlen), med_rel = median(rel))), digits = 3)
w <- wilcox.test(rel ~ sep, data = d)
cat(sprintf("Wilcoxon rel-difference, separated vs not: p = %.4f\n", w$p.value))
cat("SIGNIFICANT (separated have bigger rate gap) => long-branch attraction, not ancestry\n")
cat("NS => separation is independent of rate gap => real topological signal\n")

cat(sprintf("\nlogistic: sep ~ rel gap\n"))
print(summary(glm(sep ~ rel, data = d, family = binomial))$coefficients, digits = 3)

## does separation still hold in the LOW rate-gap half?
lo <- filter(d, rel <= median(rel))
cat(sprintf("\nlow rate-gap half only: %d/%d separated (%.2f)\n",
            sum(lo$sep), nrow(lo), mean(lo$sep)))
