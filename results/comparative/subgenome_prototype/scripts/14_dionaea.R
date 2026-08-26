#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"
k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS>0, dS<5, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
cn <- tm %>% count(nep_gene, genome, name="copies")

hr("d2/d1 — self-normalising, no speciation reference needed")
cat("  Shared A/B: both matched distances are orthologous -> ratio at the floor.\n")
cat("  Core Drosera calibrate that floor at 1.27-1.28.\n\n")
md <- k %>% filter(sp1!=sp2, sp1!=NEP, sp2!=NEP) %>%
  mutate(a=pmin(sp1,sp2), b=pmax(sp1,sp2),
         ga=ifelse(sp1==a,seq1,seq2), gb=ifelse(sp1==a,seq2,seq1)) %>%
  inner_join(cn %>% rename(a=genome, ka=copies), by=c("anchor"="nep_gene","a")) %>%
  inner_join(cn %>% rename(b=genome, kb=copies), by=c("anchor"="nep_gene","b")) %>%
  filter(ka==2, kb==2) %>%
  group_by(a,b,anchor) %>% filter(n()==4) %>%
  group_modify(function(d,key){
    ia<-sort(unique(d$ga)); ib<-sort(unique(d$gb)); M<-matrix(NA_real_,2,2)
    for(r in 1:2) for(cc in 1:2){v<-d$dS[d$ga==ia[r]&d$gb==ib[cc]]; if(length(v)) M[r,cc]<-v[1]}
    if(any(is.na(M))) return(tibble())
    p <- if (M[1,1]+M[2,2] <= M[1,2]+M[2,1]) c(M[1,1],M[2,2]) else c(M[1,2],M[2,1])
    tibble(d1=min(p), d2=max(p))
  }) %>% ungroup() %>% mutate(r = d2/d1)

res <- md %>% group_by(a,b) %>% filter(n()>=20) %>%
  summarise(loci=n(), d1=round(median(d1),3), d2=round(median(d2),3),
            ratio=round(median(d2/d1),3),
            lo=round(quantile(d2/d1,.25),2), hi=round(quantile(d2/d1,.75),2),
            has_dio = any(grepl("Dionaea", c(a,b))), .groups="drop") %>%
  arrange(ratio)
print(as.data.frame(res %>% select(-has_dio)), row.names=FALSE)

core <- c("Drosera_binata","Drosera_paradoxa","Drosera_scorpioides")
fl <- res %>% filter(a %in% core, b %in% core) %>% pull(ratio)
cat(sprintf("\n  core-Drosera floor: %s   mean %.3f\n",
            paste(fl, collapse=", "), mean(fl)))
dio <- res %>% filter(has_dio)
if (nrow(dio)) {
  cat("\n  --- Dionaea pairs ---\n")
  print(as.data.frame(dio %>% select(-has_dio)), row.names=FALSE)
  cat(sprintf("\n  Dionaea mean ratio %.3f vs floor %.3f\n",
              mean(dio$ratio), mean(fl)))
  cat("  At the floor  = Dionaea's X and Y are the SAME A and B as core Drosera.\n")
  cat("  Well above it = Dionaea's donors are a different pair. That is the answer\n")
  cat("                  to the question the whole session has been circling.\n")
  for (i in seq_len(nrow(dio))) {
    d <- md %>% filter(a==dio$a[i], b==dio$b[i])
    w <- suppressWarnings(wilcox.test(d$r, mu=mean(fl)))
    cat(sprintf("    %s / %s: n=%d, Wilcoxon vs floor p = %.3g\n",
                sub("Drosera_","",dio$a[i]), sub("Drosera_","",dio$b[i]),
                nrow(d), w$p.value))
  }
} else cat("\n  No Dionaea pair reached 20 loci with 2x2 copies — report the counts.\n")
write_csv(res, "trees/qc/d2_over_d1.csv")
