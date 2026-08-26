#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); NSIM <- 2000
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"; DIO <- "Dionaea_muscipula"

k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS>0, dS<5, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
key <- tm %>% select(anchor=nep_gene, k1=tip, genome, chr) %>% distinct() %>%
       mutate(region=sub("-.*$","",anchor))
DD <- bind_rows(k %>% transmute(anchor, ga=seq1, gb=seq2, d=dS),
                k %>% transmute(anchor, ga=seq2, gb=seq1, d=dS))
nep <- k %>% filter(xor(sp1==NEP, sp2==NEP)) %>%
       transmute(anchor, N=ifelse(sp1==NEP, seq1, seq2)) %>% distinct()

hr("0. how many 3-chromosome loci does each species have?")
tri <- key %>% filter(genome != NEP, genome != DIO) %>%
  group_by(genome, anchor, region) %>%
  summarise(nchr=n_distinct(chr), g=list(k1), c=list(chr), .groups="drop") %>%
  filter(nchr == 3, lengths(g) == 3) %>% inner_join(nep, by="anchor")
print(as.data.frame(tri %>% count(genome, name="loci")), row.names=FALSE)

hr("1. four-point within each species' own three copies")
res <- bind_rows(lapply(seq_len(nrow(tri)), function(i) {
  g <- tri$g[[i]]; cc <- tri$c[[i]]; a <- tri$anchor[i]; N <- tri$N[i]
  gp <- function(x,y){v <- DD$d[DD$anchor==a & DD$ga==x & DD$gb==y]; if(length(v)) v[1] else NA}
  m <- outer(1:3,1:3,Vectorize(function(p,q) if(p==q) NA else gp(g[p],g[q])))
  dn <- vapply(1:3, function(p) gp(g[p],N), numeric(1))
  if (any(is.na(dn)) || any(is.na(m[upper.tri(m)]))) return(NULL)
  s <- c(m[1,2]+dn[3], m[1,3]+dn[2], m[2,3]+dn[1])
  w <- which.min(s); pr <- list(c(1,2),c(1,3),c(2,3))[[w]]; og <- setdiff(1:3,pr)
  tibble(genome=tri$genome[i], region=tri$region[i], anchor=a,
         d_sis=m[pr[1],pr[2]], d_out=mean(c(m[pr[1],og], m[pr[2],og])),
         ratio=m[pr[1],pr[2]]/mean(c(m[pr[1],og], m[pr[2],og])),
         sis=paste(sort(cc[pr]), collapse="|"), support=(sort(s)[2]-sort(s)[1])/mean(s))
}))
cat(sprintf("  resolved: %d\n\n", nrow(res)))

set.seed(1)
star <- bind_rows(lapply(split(res, res$genome), function(d) {
  pool <- k$dS[k$sp1==d$genome[1] & k$sp2==d$genome[1]]
  npool <- k$dS[xor(k$sp1==NEP, k$sp2==NEP) &
                (k$sp1==d$genome[1] | k$sp2==d$genome[1])]
  if (length(pool) < 30 || length(npool) < 30) return(NULL)
  sim <- vapply(seq_len(NSIM), function(i) {
    r <- vapply(seq_len(nrow(d)), function(j) {
      mm <- sample(pool, 3, replace=TRUE); nn <- sample(npool, 3, replace=TRUE)
      s <- c(mm[1]+nn[3], mm[2]+nn[2], mm[3]+nn[1]); w <- which.min(s)
      mm[w]/mean(mm[-w])
    }, numeric(1)); median(r)
  }, numeric(1))
  tibble(genome=d$genome[1], null_ratio=round(median(sim),3),
         null_lo=round(quantile(sim,.025),3))
}))
out <- res %>% group_by(genome) %>%
  summarise(triplets=n(), d_sis=round(median(d_sis),3), d_out=round(median(d_out),3),
            ratio=round(median(ratio),3), .groups="drop") %>%
  left_join(star, by="genome") %>%
  mutate(verdict = case_when(is.na(null_ratio) ~ "n/a",
                             ratio < null_lo ~ "AAB (sisters shallower than star)",
                             TRUE ~ "consistent with ABC / star")) %>%
  arrange(ratio)
print(as.data.frame(out), row.names=FALSE)
cat("\n  The star null reproduces the selection bias: three draws from the\n")
cat("  species' own dS pool, four-point applied, same statistic computed.\n")

hr("2. is the SAME chromosome pair sister across loci? (AAB predicts yes)")
cat("  Under a star the sister pair is random -> uniform over the 3 options.\n\n")
cons <- bind_rows(lapply(split(res, paste(res$genome, res$region)), function(d) {
  if (nrow(d) < 8) return(NULL)
  tb <- sort(table(d$sis), decreasing=TRUE)
  tibble(genome=d$genome[1], region=d$region[1], n=nrow(d), n_pairs=length(tb),
         top=names(tb)[1], top_n=as.integer(tb[1]),
         frac=round(tb[1]/nrow(d),3),
         p=signif(binom.test(as.integer(tb[1]), nrow(d), 1/3)$p.value, 3))
}))
if (nrow(cons)) {
  cons$p_adj <- round(p.adjust(cons$p, "BH"), 4)
  print(as.data.frame(cons %>% mutate(genome=sub("Drosera_","D_",genome),
        region=sub("_dom","",region), top=gsub("_hap1|_collapsed","",top)) %>%
        arrange(p)), row.names=FALSE)
  cat(sprintf("\n  regions with p_adj < 0.05: %d of %d\n", sum(cons$p_adj<0.05), nrow(cons)))
} else cat("  no region reached 8 triplets\n")

hr("3. cross-check: do the sisters share a delta side?")
cat("  If AAB, the sister pair should be the two same-side segments.\n")
seg <- read_csv("trees/qc/segment_delta_calibrated.csv", show_col_types=FALSE) %>%
  mutate(side=case_when(med_delta < -0.08 ~ "A", med_delta > 0.08 ~ "B", TRUE ~ "amb"))
sp <- res %>% separate(sis, c("cA","cB"), sep="\\|", remove=FALSE) %>%
  left_join(seg %>% select(genome, region, cA=chr, sA=side), by=c("genome","region","cA")) %>%
  left_join(seg %>% select(genome, region, cB=chr, sB=side), by=c("genome","region","cB")) %>%
  filter(!is.na(sA), !is.na(sB), sA!="amb", sB!="amb")
print(as.data.frame(sp %>% group_by(genome) %>%
  summarise(n=n(), same_side=round(mean(sA==sB),3),
            ratio_same=round(median(ratio[sA==sB]),3),
            ratio_cross=round(median(ratio[sA!=sB]),3), .groups="drop")), row.names=FALSE)
cat("\n  regia gave 0.65 same-side vs 1.02 cross-side. If the core three show no\n")
cat("  such gap, their three subgenomes really are equidistant.\n")
write_csv(res, "trees/qc/triplet_all_species.csv")
