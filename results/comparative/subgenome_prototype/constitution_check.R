#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# If a species' own copies form a clade in the gene tree, they share an MRCA with each
# Dionaea copy and are FORCED to vote identically -- which alone depletes mixed
# compositions and inflates pure ones. Measure how often that happens, then redo the
# constitution spectrum on loci where the species' copies are NOT monophyletic.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ape); library(ggplot2); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
SPORD <- c("regia","binata","paradoxa","scorpioides","capensis")

ref <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE) %>%
  filter(anchor == "Nepenthes_gracilis")
v <- read_csv("tract_votes_blocks7.csv", show_col_types = FALSE) %>%
  mutate(species = sub("Drosera_","",sp))
orient <- v %>% distinct(pair, chrA, chrB) %>%
  left_join(ref %>% transmute(pair = pair_lab, winner), by = "pair") %>%
  mutate(sgA_is_X = chrA == winner)
V <- v %>% left_join(select(orient, pair, sgA_is_X), by = "pair") %>%
  mutate(g = ifelse((vote == "A") == sgA_is_X, "X", "Y"))

## per (anchor, species): are that species' own copies a clade?
meta <- read_tsv("wgd7/tip_meta.tsv", show_col_types = FALSE)
key <- function(x) gsub("@","_", gsub("['\"]","",x))
mono <- bind_rows(lapply(list.files("wgd7/tre","\\.tre$",full.names=TRUE), function(f){
  tr <- tryCatch(read.tree(f), error=function(e) NULL); if (is.null(tr)||Ntip(tr)<4) return(NULL)
  m <- meta[match(key(tr$tip.label), key(meta$tip)),]; if (any(is.na(m$genome))) return(NULL)
  nep <- tr$tip.label[m$genome=="Nepenthes_gracilis"]; if (length(nep)!=1) return(NULL)
  tr <- tryCatch(root(tr,nep,resolve.root=TRUE), error=function(e) NULL); if (is.null(tr)) return(NULL)
  m <- meta[match(key(tr$tip.label), key(meta$tip)),]
  a <- sub("\\.tre$","",basename(f))
  do.call(rbind, lapply(unique(m$genome[grepl("^Drosera_",m$genome)]), function(gg){
    tp <- tr$tip.label[m$genome == gg]; if (length(tp) < 2) return(NULL)
    data.frame(anchor=a, species=sub("Drosera_","",gg), k=length(tp),
               sisters=is.monophyletic(tr, tp), stringsAsFactors=FALSE) })) }))

cat("=== how often are a species' own copies a clade (votes forced identical)? ===\n")
print(as.data.frame(mono %>% group_by(species, k) %>%
  summarise(loci=n(), pct_clade=round(100*mean(sisters)), .groups="drop") %>%
  filter(k <= 4) %>% arrange(match(species,SPORD), k)), row.names=FALSE)

K <- V %>% group_by(species, pair, anchor) %>%
  summarise(k=n(), x=sum(g=="X"), .groups="drop") %>%
  left_join(mono %>% select(anchor, species, sisters), by=c("anchor","species")) %>%
  mutate(sisters = ifelse(is.na(sisters), FALSE, sisters))

cat("\n=== composition split by whether the copies are sisters ===\n")
print(as.data.frame(K %>% filter(k == 2) %>%
  mutate(comp = ifelse(x == 1, "mixed 1X:1Y", "pure")) %>%
  group_by(species, sisters) %>%
  summarise(loci=n(), pct_mixed=round(100*mean(comp=="mixed 1X:1Y")), .groups="drop") %>%
  arrange(match(species,SPORD), sisters)), row.names=FALSE)
cat("sisters=TRUE must show 0% mixed -- that is the artifact, by construction\n")

spec <- function(d, lab) {
  cp <- d %>% filter(k>=2, k<=4) %>% group_by(species,k) %>%
    mutate(pX=sum(x)/sum(k)) %>% group_by(species,k,x) %>%
    summarise(n=n(), pX=first(pX), .groups="drop") %>% group_by(species,k) %>%
    mutate(tot=sum(n), obs=n/tot, exp=dbinom(x,k,first(pX)),
           label=sprintf("%dX:%dY", x, k-x)) %>% ungroup()
  cat(sprintf("\n=== departure from independent sampling, %s ===\n", lab))
  print(as.data.frame(cp %>% group_by(species,k) %>%
    summarise(loci=sum(n), chisq=sum((n-tot*exp)^2/pmax(tot*exp,1e-9)),
              df=n()-2, p=ifelse(df>0, signif(pchisq(chisq,df,lower.tail=FALSE),3), NA),
              .groups="drop") %>% filter(loci >= 20) %>%
    arrange(match(species,SPORD), k)), row.names=FALSE)
  cp
}
spec(K, "ALL loci")
cpi <- spec(filter(K, !sisters), "copies NOT sisters (informative loci only)")

p1 <- cpi %>% filter(species %in% SPORD) %>% mutate(species=factor(species,SPORD)) %>%
  ggplot(aes(label, obs)) +
  geom_col(fill="#5B4EA8", width=.65) +
  geom_point(aes(y=exp), colour="#C0392B", size=1.8) +
  facet_grid(species ~ paste0(k," copies"), scales="free_x", space="free_x") +
  labs(title="Genome constitution, forced-vote loci removed",
       subtitle="only loci where the species' own copies are NOT sisters; red = independent-sampling expectation",
       x="composition (X copies : Y copies)", y="fraction of loci") +
  theme_bw(base_size=9)
ggsave("FIG37_constitution_corrected.png", p1, width=11, height=9, dpi=175, device=agg_png)
cat("\nWROTE: FIG37_constitution_corrected.png\n")
