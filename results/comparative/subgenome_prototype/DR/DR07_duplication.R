#!/usr/bin/env Rscript
# ============================================================================
# DR07 — IS THE A-DUPLICATION SHARED OR INDEPENDENT?
#
# At a locus where species X and species Y EACH have exactly two A copies,
# the four-point on those four tips answers it directly:
#
#   (X1,X2)|(Y1,Y2)  -> each species' copies are each other's closest
#                       relatives -> INDEPENDENT duplications
#   (X1,Y1)|(X2,Y2)  |
#   (X1,Y2)|(X2,Y1)  |-> copies pair ACROSS species -> SHARED duplication
#
# Copy 1 and copy 2 are arbitrary labels within a species, so the last two
# resolutions are the SAME statement. NULL IS THEREFORE 1/3 independent and
# 2/3 shared, not 1/3 each. Testing "independent" against 1/3 is the test.
#
# WHY IT MATTERS BEYOND PHASING
#   A SHARED A-duplication between regia and core Drosera is a shared derived
#   character -- evidence they are sisters, i.e. regia INSIDE Drosera.
#   INDEPENDENT duplications are consistent with regia+Dionaea.
#   So this is a second, independent line on the topology the concatenated
#   and coalescent trees disagree about.
#
# RATE-FREE: Buneman's condition uses only the ranking of three distance sums,
# and each tip enters each sum exactly once.
#
# CAPENSIS: 12-ploid with a recent WGD at dS ~0.106. Its "two A copies" may be
# recent duplicates rather than the ancestral A1/A2, which would read as
# INDEPENDENT for the wrong reason. Recent duplicates are collapsed first
# (dS < 0.25, script 15's threshold) and capensis is also reported separately.
#
# IN   DR/out/DR02_gene_delta.csv, DR/out/DR02_segments.csv,
#      DR/out/pairwise_ks.csv
# OUT  DR/out/DR07_duplication.csv
# FIG  DR/fig/DR07_1_shared.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))
COLLAPSE <- 0.25; GAPCUT <- 0.02
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

hr("1. LOCI WITH EXACTLY TWO A COPIES")
K <- read_csv("DR/out/pairwise_ks.csv", show_col_types=FALSE) %>%
     filter(!is.na(dS), dS>=0, dS<5, codons>=100)
DD <- bind_rows(K %>% transmute(locus=anchor, a=seq1, b=seq2, d=dS),
                K %>% transmute(locus=anchor, a=seq2, b=seq1, d=dS))
dmap <- setNames(DD$d, paste(DD$locus, DD$a, DD$b))
gd <- function(l,x,y) ifelse(x==y, 0, unname(dmap[paste(l,x,y)]))
G  <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE)
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE)
lab <- SS %>% filter(label %in% c("A","B")) %>%
       select(genome, region, chr, mb_lo, mb_hi, label)
GL <- G %>% mutate(mb=mid/1e6) %>%
  inner_join(lab, by=c("genome","region","chr"), relationship="many-to-many") %>%
  filter(mb>=mb_lo, mb<=mb_hi) %>% group_by(tip) %>% slice(1) %>% ungroup() %>%
  mutate(sp=sub("Drosera_","",genome)) %>% filter(label=="A")

cat("  collapsing recent duplicates at dS <", COLLAPSE, "(capensis WGD is at ~0.106)\n")
GL <- GL %>% group_by(locus, sp) %>% group_modify(function(d, key) {
  if (nrow(d) <= 1) return(d)
  n <- nrow(d); cl <- seq_len(n)
  for (i in 1:(n-1)) for (j in (i+1):n) {
    dd <- gd(key$locus, d$tip[i], d$tip[j])
    if (!is.na(dd) && dd < COLLAPSE) cl[cl==cl[j]] <- cl[i]
  }
  d[!duplicated(cl), ]
}) %>% ungroup()
two <- GL %>% count(locus, sp, name="k") %>% filter(k==2)
cat("\n  loci with exactly 2 A copies after collapse:\n")
print(as.data.frame(two %>% count(sp, name="loci")), row.names=FALSE)

hr("2. FOUR-POINT on 2 copies vs 2 copies")
cat("  null: 1/3 independent, 2/3 shared (two of three resolutions mean shared)\n\n")
SP <- sort(unique(two$sp))
res <- bind_rows(lapply(seq_along(SP), function(i)
  bind_rows(lapply(SP[-(1:i)], function(s2) {
    s1 <- SP[i]
    L <- intersect(two$locus[two$sp==s1], two$locus[two$sp==s2])
    if (!length(L)) return(NULL)
    v <- vapply(L, function(l) {
      x <- GL$tip[GL$locus==l & GL$sp==s1]; y <- GL$tip[GL$locus==l & GL$sp==s2]
      if (length(x)!=2 || length(y)!=2) return(NA_character_)
      S <- c(gd(l,x[1],x[2]) + gd(l,y[1],y[2]),
             gd(l,x[1],y[1]) + gd(l,x[2],y[2]),
             gd(l,x[1],y[2]) + gd(l,x[2],y[1]))
      if (any(is.na(S))) return(NA_character_)
      o <- sort(S)
      if ((o[2]-o[1])/mean(S) < GAPCUT) return("near-tie")
      c("independent","shared","shared")[which.min(S)]
    }, character(1))
    v <- v[!is.na(v) & v!="near-tie"]
    if (length(v) < 20) return(NULL)
    k <- sum(v=="independent"); n <- length(v)
    bt <- binom.test(k, n, 1/3)
    tibble(pair=paste(s1, s2, sep="-"), n=n,
           independent=round(k/n,3), shared=round(1-k/n,3),
           CI=sprintf("[%.2f,%.2f]", bt$conf.int[1], bt$conf.int[2]),
           p_vs_null=signif(bt$p.value,3),
           verdict=ifelse(bt$p.value>0.05, "= null",
                   ifelse(k/n > 1/3, "INDEPENDENT", "SHARED")))
  }))))
print(as.data.frame(res %>% arrange(desc(independent))), row.names=FALSE)
cat("\n  independent > 1/3 and significant -> each lineage duplicated its own A.\n")
cat("  independent < 1/3 -> copies pair across species -> SHARED duplication,\n")
cat("    which would be a shared derived character grouping those species.\n")
write_csv(res, "DR/out/DR07_duplication.csv")

hr("3. THE PAIR THAT BEARS ON THE TOPOLOGY")
rg <- res %>% filter(grepl("regia", pair))
cat("  regia vs each core species. SHARED here would group regia WITH core\n")
cat("  Drosera and argue against the regia+Dionaea topology.\n\n")
print(as.data.frame(rg), row.names=FALSE)
cc <- res %>% filter(!grepl("regia", pair), !grepl("capensis", pair))
cat("\n  core-vs-core, for comparison (do THEY share a duplication?):\n")
print(as.data.frame(cc), row.names=FALSE)

hr("4. CAPENSIS, separately")
cat("  12-ploid. Even after collapsing at dS < 0.25 its copy history is the\n")
cat("  most complex, so its rows are read with more caution.\n\n")
print(as.data.frame(res %>% filter(grepl("capensis", pair))), row.names=FALSE)

p <- ggplot(res, aes(reorder(pair, independent), independent,
                     fill=grepl("regia", pair))) +
  geom_hline(yintercept=1/3, linetype="dashed", colour="grey40") +
  geom_col(width=0.62) +
  geom_text(aes(label=sprintf("%.2f  (n=%d)", independent, n)),
            hjust=-0.12, size=3) +
  annotate("text", x=0.6, y=0.345, hjust=0, size=3, colour="grey35",
           label="1/3 = null") +
  scale_fill_manual(values=c("FALSE"="grey70","TRUE"="#D85A30"), guide="none") +
  coord_flip(ylim=c(0,1)) +
  labs(title="DR07 - is the A-duplication shared between species?",
       subtitle=paste0("Four-point on 2 A copies vs 2 A copies. Null is 1/3, ",
                       "since two of three resolutions mean 'shared'.\n",
                       "Above the line = independent duplications. ",
                       "Orange = comparisons involving regia."),
       x=NULL, y="fraction of loci resolving as INDEPENDENT duplication") +
  theme_minimal(11) + theme(plot.subtitle=element_text(size=8.5, colour="grey30"))
suppressWarnings(ggsave("DR/fig/DR07_1_shared.pdf", p, width=9, height=5))

hr("VERDICT")
cat("  ALL pairs independent -> every lineage duplicated its own A separately.\n")
cat("    Phasing A1/A2 across species is impossible, and the AAB constitution\n")
cat("    was reached more than once from the same progenitor pool.\n")
cat("  regia-core SHARED but others independent -> regia and core Drosera share\n")
cat("    a duplication, which is a shared derived character and argues regia\n")
cat("    belongs INSIDE Drosera. That would conflict with DR05/DR06.\n")
cat("  ALL shared -> one ancestral A-duplication; A1/A2 are phasable and the\n")
cat("    trees should be rebuilt with them separated.\n")
