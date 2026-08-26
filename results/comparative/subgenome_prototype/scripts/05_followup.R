#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(ape); library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd())
SET <- commandArgs(TRUE)[1]; if (is.na(SET)) SET <- "cds"
QC <- file.path(BASE,"trees","qc"); TRE <- file.path(BASE,"trees",paste0("tre_",SET))
NEP <- "Nepenthes_gracilis"
GEN <- c("Dionaea_muscipula","Drosera_binata","Drosera_capensis","Drosera_paradoxa",
         "Drosera_regia","Drosera_scorpioides","Nepenthes_gracilis")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))

trees <- read_csv(file.path(QC,paste0("tree_qc_",SET,".csv")), show_col_types=FALSE)
tips  <- read_csv(file.path(QC,paste0("tip_sides_",SET,".csv")), show_col_types=FALSE)
sp    <- read_csv(file.path(QC,paste0("species_spanning_",SET,".csv")), show_col_types=FALSE)

## fractionation partition (v2 §1.1). Prefer the file if present.
X <- c("chr1_sg1_s5","chr2_sg1_s3","chr3_sg1_s7","chr4_sg2_s12",
       "chr5_sg2_s16","chr6_sg2_s8","chr7_sg1_s1","chr8_sg1_s4")
Y <- c("chr1_sg2_s9","chr2_sg2_s11","chr3_sg2_s15","chr4_sg1_s2",
       "chr5_sg1_s10","chr6_sg1_s6","chr7_sg2_s13","chr8_sg2_s14")
pf <- file.path(BASE,"provisional_partition.csv")
if (file.exists(pf)) {
  pp <- read_csv(pf, show_col_types=FALSE); cat("provisional_partition.csv columns: ",
    paste(names(pp), collapse=", "), "\n", sep="")
}
xy_of <- function(chr) ifelse(chr %in% X, "X", ifelse(chr %in% Y, "Y", NA_character_))
dio_chr <- function(tip) sub("^Dionaea_muscipula_(chr[0-9]+_sg[0-9]+_s[0-9]+)-.*$","\\1",tip)

## ---- re-read trees for PENDANT edge lengths --------------------------------
files <- list.files(TRE, pattern="\\.treefile$", full.names=TRUE)
pend <- bind_rows(lapply(files, function(f) {
  tr <- tryCatch(read.tree(f), error=function(e) NULL); if (is.null(tr)) return(NULL)
  i <- match(seq_len(Ntip(tr)), tr$edge[,2])
  tibble(anchor=sub("\\.treefile$","",basename(f)), tip=tr$tip.label,
         pend=tr$edge.length[i])
}))
tips <- tips %>% left_join(pend, by=c("anchor","tip"))

hr("1. proper LBA test — PENDANT edge of the lone basal tip")
cat("  cumulative r2t is confounded (basal tips traverse fewer edges).\n")
cat("  Pendant edge is not.\n\n")
lone <- tips %>% group_by(anchor) %>%
  mutate(nL=sum(side=="L"), nR=sum(side=="R"),
         is_lone=(side=="L"&nL==1)|(side=="R"&nR==1),
         pend_rank=rank(-pend), n_ing=n(), pend_z=(pend-mean(pend))/sd(pend)) %>%
  ungroup()
lo <- lone %>% filter(is_lone)
cat(sprintf("  lone tip has the LONGEST pendant edge : %.1f%%  (chance %.1f%%)\n",
            100*mean(lo$pend_rank==1, na.rm=TRUE), 100*mean(1/lo$n_ing)))
cat(sprintf("  median pendant z-score of lone tip    : %+.3f  (0 = no excess)\n",
            median(lo$pend_z, na.rm=TRUE)))
cat(sprintf("  Wilcoxon lone vs nested pendant length: p = %.3g\n",
            wilcox.test(pend ~ is_lone, data=lone)$p.value))

hr("2. who owns the lone basal tip — observed vs copy-number expectation")
exp_own <- tips %>% group_by(anchor, genome) %>% summarise(k=n(), .groups="drop") %>%
  left_join(tips %>% count(anchor, name="n_ing"), by="anchor") %>%
  filter(anchor %in% (lo$anchor)) %>%
  group_by(genome) %>% summarise(expected=sum(k/n_ing), .groups="drop")
obs_own <- lo %>% count(genome, name="observed")
own <- full_join(obs_own, exp_own, by="genome") %>%
  mutate(across(c(observed,expected), ~replace_na(.,0)),
         ratio=round(observed/expected,2),
         frac_k1=sp$genome %>% {sapply(genome, function(g)
           round(mean(sp$k[sp$genome==g]==1),3))}) %>% arrange(desc(ratio))
print(as.data.frame(own), row.names=FALSE)
cat("\n  Under shared A/B, only species that reliably keep BOTH copies can be\n")
cat("  lone-basal: if every single-copy Drosera kept the same lineage, the\n")
cat("  other Dionaea copy is alone by necessity. Predicts ratio > 1 exactly\n")
cat("  for the low-frac_k1 species.\n")

hr("3. is the basal Dionaea copy consistently X or Y?")
dl <- lo %>% filter(genome=="Dionaea_muscipula") %>%
  mutate(chr=dio_chr(tip), xy=xy_of(chr))
cat(sprintf("  lone-basal Dionaea tips: %d   assigned to partition: %d\n",
            nrow(dl), sum(!is.na(dl$xy))))
print(table(dl$xy, useNA="ifany"))
if (sum(!is.na(dl$xy)) > 20) {
  b <- binom.test(sum(dl$xy=="X", na.rm=TRUE), sum(!is.na(dl$xy)))
  cat(sprintf("\n  binomial vs 50/50: p = %.3g\n", b$p.value))
  cat("  ~50/50 = consistent with the fractionation story (which lineage the\n")
  cat("           Drosera happened to keep varies by locus).\n")
  cat("  strongly skewed = one Dionaea subgenome is systematically basal, a\n")
  cat("           much stronger and more surprising claim. Follow up either way.\n")
}
cat("\n  per Dionaea chromosome pair:\n")
print(as.data.frame(dl %>% count(chr, xy) %>% arrange(chr)), row.names=FALSE)

hr("4. rate control, DE-CONFOUNDED for species")
rt <- tips %>% group_by(anchor, genome) %>%
  summarise(rate_gap=if (n()>=2) max(r2t)-min(r2t) else NA_real_,
            pend_gap=if (n()>=2) max(pend)-min(pend) else NA_real_,
            mean_r2t=mean(r2t), .groups="drop")
d <- sp %>% inner_join(rt, by=c("anchor","genome")) %>%
  left_join(trees %>% select(anchor,sup_min,balance,n_species), by="anchor") %>%
  filter(genome != NEP, !is.na(spans), !is.na(rate_gap))
cat("  a) pooled WITH species as a covariate:\n")
print(round(summary(glm(spans ~ scale(rate_gap)+scale(mean_r2t)+genome,
      data=d, family=binomial))$coefficients[1:3,], 4))
cat("\n  b) within each species separately (the honest version):\n")
per <- lapply(split(d, d$genome), function(x) {
  if (nrow(x) < 40) return(NULL)
  g <- glm(spans ~ scale(rate_gap)+scale(mean_r2t), data=x, family=binomial)
  cf <- summary(g)$coefficients
  tibble(genome=x$genome[1], n=nrow(x),
         b_rategap=round(cf["scale(rate_gap)","Estimate"],3),
         p_rategap=signif(cf["scale(rate_gap)","Pr(>|z|)"],3),
         b_meanr2t=round(cf["scale(mean_r2t)","Estimate"],3),
         p_meanr2t=signif(cf["scale(mean_r2t)","Pr(>|z|)"],3))
})
print(as.data.frame(bind_rows(per)), row.names=FALSE)
cat("\n  Positive b_rategap WITHIN a species = LBA for that species.\n")

hr("5. what the clean-subset filter selects for")
cl <- trees %>% mutate(clean = sup_min>=80 & balance>=0.3 & n_species==6)
cat(sprintf("  clean: %d   rest: %d\n\n", sum(cl$clean), sum(!cl$clean)))
print(as.data.frame(cl %>% group_by(clean) %>%
  summarise(n=n(), med_ntips=median(n_ingroup), med_sup=median(sup_min,na.rm=TRUE),
            med_balance=round(median(balance),3),
            med_nepratio=round(median(nep_ratio,na.rm=TRUE),2), .groups="drop")),
  row.names=FALSE)
cat("\n  If clean trees are systematically bigger, the 5/6-species result may\n")
cat("  just be the trees with most information. Check n before celebrating.\n")
