#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE)
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"; DIO <- "Dionaea_muscipula"; REG <- "Drosera_regia"
X <- c("chr1_sg1_s5","chr2_sg1_s3","chr3_sg1_s7","chr4_sg2_s12",
       "chr5_sg2_s16","chr6_sg2_s8","chr7_sg1_s1","chr8_sg1_s4")

k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS>0, dS<5, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)

hr("0. join diagnostic — which meta column matches the ks seq names?")
sq <- unique(c(k$seq1, k$seq2))
for (nm in c("tip","gene")) {
  cat(sprintf("  %-22s %.4f\n", nm, mean(sq %in% tm[[nm]])))
  cat(sprintf("  %-22s %.4f\n", paste0(nm," @->_"), mean(sq %in% gsub("@","_",tm[[nm]]))))
}
cands <- list(tip=tm$tip, tip_=gsub("@","_",tm$tip), gene=tm$gene, gene_=gsub("@","_",tm$gene))
best <- names(which.max(vapply(cands, function(v) mean(sq %in% v), numeric(1))))
cat(sprintf("\n  using: %s\n", best))
key <- tm %>% mutate(k1 = cands[[best]]) %>% select(anchor=nep_gene, k1, genome, chr) %>% distinct()
if (max(vapply(cands, function(v) mean(sq %in% v), numeric(1))) < 0.9)
  stop("no column joins — inspect tip_meta and pairwise_ks headers")

WI <- k %>% filter(sp1==sp2, sp1!=NEP) %>% transmute(anchor, sp=sp1, g1=seq1, g2=seq2, dS) %>%
  left_join(key %>% rename(g1=k1, chr1=chr), by=c("anchor","g1","sp"="genome")) %>%
  left_join(key %>% rename(g2=k1, chr2=chr), by=c("anchor","g2","sp"="genome")) %>%
  mutate(same_chr = chr1 == chr2)

hr("1. is the recent mode TANDEM? within-species pairs by same/different chromosome")
print(as.data.frame(WI %>% filter(!is.na(same_chr)) %>% group_by(sp, same_chr) %>%
  summarise(pairs=n(), p25=round(quantile(dS,.25),3), med=round(median(dS),3),
            p75=round(quantile(dS,.75),3), .groups="drop")), row.names=FALSE)
cat("\n  If regia's same_chr pairs sit at ~0.23 and diff_chr much higher, the\n")
cat("  'recent WGD' is array duplication and can simply be dropped.\n")

hr("2. regia copy number, before and after dropping same-chromosome copies")
rc <- key %>% filter(genome==REG) %>% count(anchor, name="copies")
rc2 <- key %>% filter(genome==REG) %>% group_by(anchor) %>%
       summarise(chrs = n_distinct(chr), .groups="drop")
cat("  raw copies per locus:\n"); print(table(rc$copies))
cat("\n  distinct chromosomes per locus:\n"); print(table(rc2$chrs))
cat(sprintf("\n  loci with 3+ distinct regia chromosomes: %d\n", sum(rc2$chrs>=3)))

hr("3. per ancestral region — which regia chromosomes carry it")
rr <- key %>% filter(genome==REG) %>% mutate(region=sub("-.*$","",anchor)) %>%
      group_by(region, chr) %>% summarise(loci=n(), .groups="drop") %>%
      group_by(region) %>% slice_max(loci, n=6) %>% ungroup()
print(as.data.frame(rr %>% group_by(region) %>%
  summarise(n_chr=n(), top=paste(sub("_collapsed","",chr), loci, sep=":", collapse=" "),
            .groups="drop")), row.names=FALSE)
cat("\n  ~3 well-represented chromosomes per region = the three subgenomes.\n")

hr("4. triplet resolution within regia (four-point, Nepenthes outgroup)")
cat("  At 3-chromosome loci: which two regia copies are closest?\n\n")
DD <- bind_rows(k %>% transmute(anchor, ga=seq1, gb=seq2, d=dS),
                k %>% transmute(anchor, ga=seq2, gb=seq1, d=dS))
nep <- k %>% filter(xor(sp1==NEP, sp2==NEP)) %>%
       transmute(anchor, N=ifelse(sp1==NEP, seq1, seq2)) %>% distinct()
tri <- key %>% filter(genome==REG) %>% group_by(anchor) %>%
  filter(n_distinct(chr)>=3) %>% slice_head(n=3) %>%
  summarise(g=list(k1), c=list(chr), .groups="drop") %>%
  inner_join(nep, by="anchor")
res <- bind_rows(lapply(seq_len(nrow(tri)), function(i) {
  g <- tri$g[[i]]; cc <- tri$c[[i]]; N <- tri$N[i]; a <- tri$anchor[i]
  gp <- function(x,y) { v <- DD$d[DD$anchor==a & DD$ga==x & DD$gb==y]; if (length(v)) v[1] else NA }
  m <- outer(1:3, 1:3, Vectorize(function(p,q) if (p==q) NA else gp(g[p], g[q])))
  dn <- vapply(1:3, function(p) gp(g[p], N), numeric(1))
  if (any(is.na(dn)) || any(is.na(m[upper.tri(m)]))) return(NULL)
  s <- c(m[1,2]+dn[3], m[1,3]+dn[2], m[2,3]+dn[1])
  w <- which.min(s); pr <- list(c(1,2), c(1,3), c(2,3))[[w]]
  tibble(anchor=a, pair=paste(sort(sub("_collapsed","",cc[pr])), collapse="|"),
         out=sub("_collapsed","",cc[setdiff(1:3,pr)]),
         support=(sort(s)[2]-sort(s)[1])/mean(s))
}))
cat(sprintf("  resolved triplets: %d\n\n", nrow(res)))
if (nrow(res)) {
  print(as.data.frame(res %>% count(pair, sort=TRUE) %>% head(15)), row.names=FALSE)
  cat(sprintf("\n  median support: %.3f   unresolved (<0.05): %.2f\n",
              median(res$support), mean(res$support < 0.05)))
  write_csv(res, "trees/qc/regia_triplets.csv")
}

hr("5. each regia chromosome placed on the Dionaea A/B axis")
cat("  delta = [d(R, Dio_X) - d(R, Dio_Y)] / d(X,Y). Negative = allies with X.\n")
cat("  Two groups at 2:1 would mean regia is AAB rather than ABC.\n\n")
dio <- key %>% filter(genome==DIO) %>% mutate(xy=ifelse(chr %in% X,"X","Y")) %>%
       select(anchor, k1, xy) %>%
       pivot_wider(names_from=xy, values_from=k1, values_fn=list) %>%
       filter(lengths(X)==1, lengths(Y)==1) %>%
       mutate(DX=unlist(X), DY=unlist(Y)) %>% select(anchor, DX, DY)
dd <- DD %>% rename(DX=ga, DY=gb, dXY=d)
D <- key %>% filter(genome==REG) %>% inner_join(dio, by="anchor") %>%
  inner_join(dd, by=c("anchor","DX","DY")) %>%
  left_join(DD %>% rename(k1=ga, DX=gb, dRX=d), by=c("anchor","k1","DX")) %>%
  left_join(DD %>% rename(k1=ga, DY=gb, dRY=d), by=c("anchor","k1","DY")) %>%
  filter(!is.na(dRX), !is.na(dRY), dXY>0.05) %>%
  mutate(delta=(dRX-dRY)/dXY)
print(as.data.frame(D %>% group_by(chr) %>% filter(n()>=25) %>%
  summarise(copies=n(), med_delta=round(median(delta),3),
            fracX=round(mean(delta < 0),3), .groups="drop") %>%
  arrange(med_delta)), row.names=FALSE)
cat("\n  Compare to v2 11.3: chr2 and chr11 were the Y-set there.\n")
write_csv(D, "trees/qc/regia_delta_by_chr.csv")

hr("6. AA'B or AAB? depth ratio within each resolved triplet")
cat("  ratio = d(sister pair) / mean d(sister, outgroup).\n")
cat("  ~0.3 = the sisters are a RECENT duplicate  -> AA'B\n")
cat("  ~1.0 = all three lineages comparably deep  -> AAB\n\n")
dep <- bind_rows(lapply(seq_len(nrow(tri)), function(i) {
  g <- tri$g[[i]]; cc <- tri$c[[i]]; N <- tri$N[i]; a <- tri$anchor[i]
  gp <- function(x,y){v <- DD$d[DD$anchor==a & DD$ga==x & DD$gb==y]; if(length(v)) v[1] else NA}
  m <- outer(1:3,1:3,Vectorize(function(p,q) if(p==q) NA else gp(g[p],g[q])))
  dn <- vapply(1:3, function(p) gp(g[p],N), numeric(1))
  if (any(is.na(dn)) || any(is.na(m[upper.tri(m)]))) return(NULL)
  s <- c(m[1,2]+dn[3], m[1,3]+dn[2], m[2,3]+dn[1])
  w <- which.min(s); pr <- list(c(1,2),c(1,3),c(2,3))[[w]]; og <- setdiff(1:3,pr)
  tibble(anchor=a, d_sis=m[pr[1],pr[2]],
         d_out=mean(c(m[pr[1],og], m[pr[2],og])),
         ratio=m[pr[1],pr[2]]/mean(c(m[pr[1],og], m[pr[2],og])),
         sis=paste(sort(sub("_collapsed","",cc[pr])),collapse="|"))
}))
cat(sprintf("  triplets: %d\n", nrow(dep)))
cat(sprintf("  median d(sisters) = %.3f   median d(sister,outgroup) = %.3f\n",
            median(dep$d_sis), median(dep$d_out)))
cat(sprintf("  median ratio = %.3f   [IQR %.3f - %.3f]\n",
            median(dep$ratio), quantile(dep$ratio,.25), quantile(dep$ratio,.75)))
cat(sprintf("  triplets with ratio < 0.5: %.2f\n", mean(dep$ratio < 0.5)))
cat("\n  by sister pair (>=5 triplets):\n")
print(as.data.frame(dep %>% group_by(sis) %>% filter(n()>=5) %>%
  summarise(n=n(), d_sis=round(median(d_sis),3), d_out=round(median(d_out),3),
            ratio=round(median(ratio),3), .groups="drop") %>% arrange(ratio)),
  row.names=FALSE)
cat("\n  The ratio is within-locus, so gene rate cancels; both terms are\n")
cat("  within-regia, so lineage rate cancels too.\n")
write_csv(dep, "trees/qc/regia_triplet_depths.csv")
