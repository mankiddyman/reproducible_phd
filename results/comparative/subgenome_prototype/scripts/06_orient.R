#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd())
SET <- commandArgs(TRUE)[1]; if (is.na(SET)) SET <- "cds"
QC <- file.path(BASE,"trees","qc")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))

trees <- read_csv(file.path(QC,paste0("tree_qc_",SET,".csv")), show_col_types=FALSE)
tips  <- read_csv(file.path(QC,paste0("tip_sides_",SET,".csv")), show_col_types=FALSE)
meta  <- read_tsv(file.path(BASE,"wgd7","tip_meta.tsv"), show_col_types=FALSE) %>%
         mutate(tip = gsub("@","_",tip)) %>% select(anchor=nep_gene, tip, chr)
tips <- tips %>% left_join(meta, by=c("anchor","tip")) %>%
        mutate(region = sub("-.*$","",anchor))

X <- c("chr1_sg1_s5","chr2_sg1_s3","chr3_sg1_s7","chr4_sg2_s12",
       "chr5_sg2_s16","chr6_sg2_s8","chr7_sg1_s1","chr8_sg1_s4")
Y <- c("chr1_sg2_s9","chr2_sg2_s11","chr3_sg2_s15","chr4_sg1_s2",
       "chr5_sg1_s10","chr6_sg1_s6","chr7_sg2_s13","chr8_sg2_s14")

## orientation: keep trees where Dionaea X and Y are on OPPOSITE sides
dio <- tips %>% filter(genome=="Dionaea_muscipula") %>%
  mutate(xy = ifelse(chr %in% X, "X", ifelse(chr %in% Y, "Y", NA))) %>%
  filter(!is.na(xy)) %>% select(anchor, side, xy) %>% distinct()
ori <- dio %>% group_by(anchor) %>%
  filter(n()==2, n_distinct(side)==2, n_distinct(xy)==2) %>%
  summarise(x_side = side[xy=="X"], .groups="drop")
cat(sprintf("trees with Dionaea X and Y on opposite sides: %d of %d (%.1f%%)\n",
            nrow(ori), nrow(trees), 100*nrow(ori)/nrow(trees)))

## every non-Dionaea tip votes X-side or Y-side
v <- tips %>% inner_join(ori, by="anchor") %>%
  filter(genome != "Dionaea_muscipula", genome != "Nepenthes_gracilis") %>%
  mutate(votes_X = side == x_side)

## analytic baseline: per tree, P(a random non-Dio tip lands X-side)
base <- v %>% group_by(anchor) %>%
  summarise(p_x = mean(votes_X), n_t = n(), .groups="drop")
v <- v %>% left_join(base, by="anchor")

zt <- function(d) {
  n <- nrow(d); if (n < 15) return(tibble(n=n, obs=NA_real_, exp=NA_real_, z=NA_real_))
  o <- sum(d$votes_X); e <- sum(d$p_x); vv <- sum(d$p_x*(1-d$p_x))
  tibble(n=n, obs=round(o/n,3), exp=round(e/n,3),
         z=round(if (vv>0) (o-e)/sqrt(vv) else NA_real_, 2))
}

hr("1. species-level skew toward Dionaea X or Y")
cat("  baseline is NOT 0.5 — it is set per tree by how many tips sit each side.\n\n")
print(as.data.frame(v %>% group_by(genome) %>% group_modify(~zt(.x)) %>% ungroup() %>%
        arrange(desc(z))), row.names=FALSE)

hr("2. same, restricted to well-supported trees (sup_min >= 80)")
keep <- trees$anchor[!is.na(trees$sup_min) & trees$sup_min >= 80]
cat(sprintf("  oriented trees kept: %d\n\n", sum(ori$anchor %in% keep)))
print(as.data.frame(v %>% filter(anchor %in% keep) %>% group_by(genome) %>%
        group_modify(~zt(.x)) %>% ungroup() %>% arrange(desc(z))), row.names=FALSE)

hr("3. support sweep — does the signal strengthen with resolution?")
print(as.data.frame(bind_rows(lapply(c(0,50,70,80,90,95), function(th) {
  k <- trees$anchor[!is.na(trees$sup_min) & trees$sup_min >= th]
  d <- v %>% filter(anchor %in% k)
  zt(d) %>% mutate(sup_min_ge = th, trees = n_distinct(d$anchor)) %>%
    select(sup_min_ge, trees, n, obs, exp, z)
}))), row.names=FALSE)
cat("\n  z rising with threshold = real signal sharpened by better trees.\n")
cat("  z flat or falling = the signal lives in the unresolved trees. Bad sign.\n")

hr("4. per Drosera chromosome — the co-occurrence result")
cat("  Consistent X-side and Y-side sets = one shared A/B partition.\n\n")
ch <- v %>% group_by(genome, region, chr) %>% group_modify(~zt(.x)) %>% ungroup() %>%
      filter(!is.na(z), n >= 20) %>% arrange(z)
cat("  --- strongest Y-side ---\n")
print(as.data.frame(head(ch, 12)), row.names=FALSE)
cat("\n  --- strongest X-side ---\n")
print(as.data.frame(tail(ch, 12)), row.names=FALSE)
cat(sprintf("\n  chromosomes tested: %d   |z| > 2: %d   |z| > 3: %d\n",
            nrow(ch), sum(abs(ch$z)>2), sum(abs(ch$z)>3)))

hr("5. LBA control on the orientation")
cat("  X is the FASTER subgenome (v2 3.2, 7/8 pairs). If this is LBA, fast\n")
cat("  Drosera tips should be pulled to the X side. Test that directly.\n\n")
g <- glm(votes_X ~ scale(r2t) + genome, data=v, family=binomial)
print(round(summary(g)$coefficients[1:2,,drop=FALSE], 4))
cat("\n  Positive scale(r2t) = fast tips go X-side = LBA. Null or negative = not.\n")

write_csv(ch, file.path(QC, paste0("chrom_xy_votes_", SET, ".csv")))
write_csv(v %>% select(anchor, region, tip, genome, chr, votes_X, p_x, r2t),
          file.path(QC, paste0("xy_votes_", SET, ".csv")))
cat(sprintf("\nWROTE: %s/{chrom_xy_votes,xy_votes}_%s.csv\n", QC, SET))

hr("6. THE key split — forced vs informative trees")
lone_anch <- tips %>% filter(genome=="Dionaea_muscipula") %>%
  group_by(anchor) %>% mutate(nL=sum(side=="L"), nR=sum(side=="R")) %>%
  filter((side=="L"&nL==1&nR==0)|(side=="R"&nR==1&nL==0)) %>% pull(anchor) %>% unique()
vv <- v %>% mutate(forced = anchor %in% lone_anch)
cat(sprintf("  forced (a Dionaea copy alone): %d trees\n",
            n_distinct(vv$anchor[vv$forced])))
cat(sprintf("  informative:                   %d trees\n\n",
            n_distinct(vv$anchor[!vv$forced])))
print(as.data.frame(vv %>% group_by(forced, genome) %>% group_modify(~zt(.x)) %>%
        ungroup() %>% arrange(forced, desc(z))), row.names=FALSE)
cat("\n  In INFORMATIVE trees: do Drosera tips ever land Y-side?\n")
print(as.data.frame(vv %>% filter(!forced) %>% group_by(anchor) %>%
        summarise(n_Y = sum(!votes_X), .groups="drop") %>% count(n_Y)), row.names=FALSE)
cat("\n  Mass at n_Y = 0 across the board = Y is a lineage no Drosera has.\n")
cat("  A spread = Drosera really are A/B mixtures. This is the question.\n")
