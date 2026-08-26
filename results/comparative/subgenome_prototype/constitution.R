#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# 1. Is forced Drosera monophyly explained by LOW COPY NUMBER (subgenome loss)?
# 2. Genome constitution: at loci where a species has exactly k copies, how do those
#    k copies split X vs Y? A fixed constitution (say AAB) gives 2:1 far more often
#    than independent per-copy sampling at the same marginal rate would.
# Orientation comes from fractionation; a global flip would swap AAB<->ABB, so read
# DIFFERENCES between species as robust and the letters themselves as conditional.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ape); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
set.seed(1)
SPORD <- c("regia","binata","paradoxa","scorpioides","capensis")

## ---- orientation: local sgA -> global X ------------------------------------
ref <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE) %>%
  filter(anchor == "Nepenthes_gracilis")
v <- read_csv("tract_votes_blocks7.csv", show_col_types = FALSE) %>%
  mutate(species = sub("Drosera_","",sp))
orient <- v %>% distinct(pair, chrA, chrB) %>%
  left_join(ref %>% transmute(pair = pair_lab, winner), by = "pair") %>%
  mutate(sgA_is_X = chrA == winner)
stopifnot(!any(is.na(orient$sgA_is_X)))
V <- v %>% left_join(select(orient, pair, sgA_is_X), by = "pair") %>%
  mutate(g = ifelse((vote == "A") == sgA_is_X, "X", "Y"))
cat(sprintf("votes %d | globally oriented | %.0f%% X overall\n", nrow(V), 100*mean(V$g=="X")))

## ---- 1. forced monophyly vs Drosera copy number ----------------------------
meta <- read_tsv("wgd7/tip_meta.tsv", show_col_types = FALSE)
key <- function(x) gsub("@","_", gsub("['\"]","",x))
fx <- bind_rows(lapply(list.files("wgd7/tre","\\.tre$",full.names=TRUE), function(f){
  tr <- tryCatch(read.tree(f), error=function(e) NULL); if (is.null(tr)||Ntip(tr)<4) return(NULL)
  m <- meta[match(key(tr$tip.label), key(meta$tip)),]; if (any(is.na(m$genome))) return(NULL)
  nep <- tr$tip.label[m$genome=="Nepenthes_gracilis"]; if (length(nep)!=1) return(NULL)
  tr <- tryCatch(root(tr,nep,resolve.root=TRUE), error=function(e) NULL); if (is.null(tr)) return(NULL)
  m <- meta[match(key(tr$tip.label), key(meta$tip)),]
  dio <- tr$tip.label[m$genome=="Dionaea_muscipula"]
  dro <- tr$tip.label[grepl("^Drosera_", m$genome)]
  if (length(dio)!=2 || length(dro)<2) return(NULL)
  tibble(anchor=sub("\\.tre$","",basename(f)), n_dro=length(dro),
         n_sp=n_distinct(m$genome[grepl("^Drosera_",m$genome)]),
         forced=is.monophyletic(tr,dro)) }))
cat("\n=== is forced monophyly driven by low Drosera copy number? ===\n")
print(as.data.frame(fx %>% mutate(bin = cut(n_dro, c(1,2,3,4,6,100),
        labels=c("2","3","4","5-6","7+"))) %>%
  group_by(bin) %>% summarise(trees=n(), pct_forced=round(100*mean(forced)), .groups="drop")))
cat("falling with copy number => subgenome LOSS explains it (your reading)\n")
cat("flat across copy number  => gene conversion is the better explanation\n")

## ---- 2. constitution spectrum ----------------------------------------------
K <- V %>% group_by(species, pair, anchor) %>%
  summarise(k = n(), x = sum(g == "X"), .groups = "drop")
cat("\n=== loci per species by observed copy number ===\n")
print(as.data.frame(K %>% count(species, k) %>%
  pivot_wider(names_from=k, values_from=n, values_fill=0) %>% arrange(match(species,SPORD))))

comp <- K %>% filter(k >= 2, k <= 4) %>%
  group_by(species, k) %>%
  mutate(pX = sum(x)/sum(k)) %>%           # that species' marginal X rate at that k
  group_by(species, k, x) %>%
  summarise(n = n(), pX = first(pX), .groups = "drop") %>%
  group_by(species, k) %>%
  mutate(tot = sum(n), obs = n/tot,
         exp = dbinom(x, k, first(pX)),
         label = sprintf("%dX:%dY", x, k - x)) %>% ungroup()
cat("\n=== observed vs independent-sampling expectation ===\n")
for (kk in 2:3) {
  cat(sprintf("\n-- %d copies --\n", kk))
  print(as.data.frame(comp %>% filter(k == kk) %>%
    transmute(species, label, loci = n, obs = round(obs,3), exp = round(exp,3),
              enrich = round(obs/exp,2)) %>% arrange(match(species,SPORD), label)),
    row.names = FALSE)
}

## chi-square: does the composition depart from independent sampling?
cat("\n=== departure from independent per-copy sampling ===\n")
print(as.data.frame(comp %>% group_by(species, k) %>%
  summarise(loci = sum(n),
            chisq = sum((n - tot*exp)^2 / pmax(tot*exp, 1e-9)),
            df = n() - 2,
            p = ifelse(df > 0, signif(pchisq(chisq, df, lower.tail=FALSE), 3), NA),
            .groups="drop") %>% arrange(match(species,SPORD), k)), row.names = FALSE)
cat("small p = copies are NOT independent => a structured constitution (AAB/ABB)\n")
cat("large p = copies look like independent draws => no fixed constitution detectable\n")

## ---- 3. modal constitution --------------------------------------------------
cat("\n=== modal composition per species per copy number ===\n")
print(as.data.frame(comp %>% group_by(species, k) %>%
  slice_max(obs, n = 1, with_ties = FALSE) %>%
  transmute(species, k, modal = label, frac_of_loci = round(obs,3), loci = n) %>%
  arrange(match(species,SPORD), k)), row.names = FALSE)

cat("\n=== global X fraction per species (regia vs the rest) ===\n")
print(as.data.frame(V %>% group_by(species) %>%
  summarise(votes=n(), fracX=round(mean(g=="X"),3),
            p=signif(binom.test(sum(g=="X"), n(), .5)$p.value,3), .groups="drop") %>%
  arrange(match(species,SPORD))), row.names = FALSE)

p1 <- comp %>% mutate(species = factor(species, SPORD)) %>%
  ggplot(aes(label, obs, fill = "observed")) +
  geom_col(width=.65) +
  geom_point(aes(y = exp, colour = "independent-sampling expectation"), size=1.8) +
  facet_grid(species ~ paste0(k, " copies"), scales="free_x", space="free_x") +
  scale_fill_manual(values=c(observed="#5B4EA8"), name=NULL) +
  scale_colour_manual(values=c(`independent-sampling expectation`="#C0392B"), name=NULL) +
  labs(title="Genome constitution: how a species' copies of one gene split between the two Dionaea subgenomes",
       subtitle="bars = observed fraction of loci; red points = expected if each copy were drawn independently at that species' marginal X rate",
       x="composition (X copies : Y copies)", y="fraction of loci") +
  theme_bw(base_size=9) + theme(legend.position="top")
ggsave("FIG36_constitution.png", p1, width=11, height=9, dpi=175, device=agg_png)
write_csv(comp, "constitution_spectrum.csv")
cat("\nWROTE: FIG36_constitution.png constitution_spectrum.csv\n")
