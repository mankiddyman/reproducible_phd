#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Does the MINIMUM cross-species dS equal the speciation distance?
#   shared A/B   -> A_X vs A_Y is an ortholog pair -> min == speciation
#   independent  -> the donor sublineages had already split -> min > speciation
# Rates cancel exactly: both quantities involve the same two species.
# Also reports MAX within-species dS per locus, which estimates the DEEPEST duplication
# and is robust to extra recent WGDs (the capensis / regia problem).
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ragg); library(patchwork)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
SPORD <- c("Dionaea_muscipula","Drosera_regia","Drosera_binata",
           "Drosera_paradoxa","Drosera_scorpioides","Drosera_capensis")
NEP <- "Nepenthes_gracilis"

k <- read_csv("ks/pairwise_ks.csv", show_col_types = FALSE) %>%
  filter(!is.na(dS), dS > 0, dS < 5, codons >= 100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types = FALSE)
cn <- tm %>% count(nep_gene, genome, name = "copies")
sc <- cn %>% filter(copies == 1)
md <- function(x, to = 3) { if (length(x) < 25) return(NA_real_)
  d <- density(x, from = 0, to = to, n = 2048); d$x[which.max(d$y)] }

## ---- reference: speciation dS from single-copy orthologs -------------------
SPEC <- k %>% filter(sp1 != sp2, sp1 != NEP, sp2 != NEP) %>%
  semi_join(sc, by = c("anchor"="nep_gene","sp1"="genome")) %>%
  semi_join(sc, by = c("anchor"="nep_gene","sp2"="genome")) %>%
  mutate(a = pmin(sp1,sp2), b = pmax(sp1,sp2)) %>%
  group_by(a,b) %>% summarise(n_sc = n(), spec = median(dS), .groups="drop")
cat("=== speciation reference (single-copy orthologs) ===\n")
print(as.data.frame(SPEC %>% transmute(a,b,n_sc,spec=round(spec,3))), row.names=FALSE)

## ---- deepest duplication per locus (robust to extra recent WGDs) -----------
DEEP <- k %>% filter(sp1 == sp2, sp1 != NEP) %>% rename(species = sp1) %>%
  group_by(anchor, species) %>%
  summarise(pairs = n(), dup_min = min(dS), dup_max = max(dS), .groups="drop")
cat("\n=== deepest within-species duplication (max over that locus's pairs) ===\n")
cat("min = most recent event, max = deepest. A species with several nested WGDs\n")
cat("shows a large gap between them; a median would sit between and mean nothing.\n\n")
print(as.data.frame(DEEP %>% group_by(species) %>%
  summarise(loci = n(), mode_min = round(md(dup_min),3), mode_max = round(md(dup_max),3),
            med_min = round(median(dup_min),3), med_max = round(median(dup_max),3),
            .groups="drop") %>% arrange(match(species,SPORD))), row.names=FALSE)

## ---- THE TEST: min cross-species dS vs speciation --------------------------
multi <- cn %>% filter(copies >= 2) %>% transmute(anchor = nep_gene, species = genome)
X <- k %>% filter(sp1 != sp2, sp1 != NEP, sp2 != NEP) %>%
  semi_join(multi, by = c("anchor","sp1"="species")) %>%
  semi_join(multi, by = c("anchor","sp2"="species")) %>%
  mutate(a = pmin(sp1,sp2), b = pmax(sp1,sp2)) %>%
  group_by(anchor, a, b) %>%
  summarise(n_cross = n(), cross_min = min(dS), cross_max = max(dS), .groups="drop") %>%
  filter(n_cross >= 4) %>%                     # need >=2 copies in each species
  left_join(SPEC, by = c("a","b")) %>% filter(!is.na(spec)) %>%
  mutate(ratio = cross_min / spec)
cat(sprintf("\nloci with >=2 copies in BOTH species: %d across %d species pairs\n",
            nrow(X), n_distinct(paste(X$a,X$b))))

RES <- X %>% group_by(a,b) %>%
  summarise(loci = n(), spec = round(first(spec),3),
            cross_min_med = round(median(cross_min),3),
            ratio_med = round(median(ratio),2),
            frac_near = round(mean(ratio > 0.75 & ratio < 1.33),3),
            .groups="drop") %>% arrange(ratio_med)
cat("\n=== min cross-species dS vs speciation ===\n")
cat("ratio ~1  => the closest cross-species copies ARE orthologs => SHARED A/B\n")
cat("ratio >1  => already separated before the merger => INDEPENDENT hybridisations\n\n")
print(as.data.frame(RES), row.names=FALSE)
cat(sprintf("\npairs with median ratio in 0.75-1.33: %d of %d\n",
            sum(RES$ratio_med > 0.75 & RES$ratio_med < 1.33), nrow(RES)))
write_csv(RES, "shared_ab_test.csv")

## ---- plots -----------------------------------------------------------------
p1 <- X %>% mutate(pair = paste(sub("Drosera_","",sub("Dionaea_","",a)), "/",
                                sub("Drosera_","",sub("Dionaea_","",b)))) %>%
  ggplot(aes(ratio)) +
  geom_histogram(binwidth = .1, boundary = 0, fill = "#5B4EA8", colour = NA) +
  geom_vline(xintercept = 1, colour = "#C0392B", linetype = 2) +
  facet_wrap(~ pair, scales = "free_y", ncol = 3) +
  coord_cartesian(xlim = c(0, 3)) +
  labs(title = "Is the closest cross-species copy pair an ortholog?",
       subtitle = "ratio = min(cross-species dS) / speciation dS.  1 (red) = shared A/B; >1 = the donor lineages had already split.",
       x = "min cross-species dS / speciation dS", y = "loci") +
  theme_bw(base_size = 8)
ggsave("FIG45_shared_ab_test.png", p1, width = 12, height = 8, dpi = 180, device = agg_png)

p2 <- DEEP %>% pivot_longer(c(dup_min, dup_max), names_to = "which", values_to = "dS") %>%
  mutate(species = factor(species, SPORD),
         which = recode(which, dup_min = "most recent event", dup_max = "deepest event")) %>%
  ggplot(aes(dS, fill = which)) +
  geom_density(alpha = .5, colour = NA, bw = .05) +
  facet_wrap(~ species, ncol = 3, scales = "free_y") +
  coord_cartesian(xlim = c(0, 2)) +
  scale_fill_manual(values = c(`most recent event` = "#d95f02", `deepest event` = "#1b9e77")) +
  labs(title = "Separating nested WGDs within a species",
       subtitle = "per locus: the smallest and largest dS among that species' own copies. Two separated peaks = two events.",
       x = "dS", y = NULL, fill = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "top")
ggsave("FIG46_nested_wgds.png", p2, width = 11, height = 7, dpi = 180, device = agg_png)
cat("\nWROTE: FIG45_shared_ab_test.png FIG46_nested_wgds.png shared_ab_test.csv\n")
