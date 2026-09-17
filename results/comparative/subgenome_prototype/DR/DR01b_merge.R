#!/usr/bin/env Rscript
# ============================================================================
# DR01b — CAN MERGING BLOCKS RESCUE THE VOTE TEST?
#
# The vote test failed because GENESPACE splits collinear tracts into many
# small blocks: n_units per region reached 10-41 against an assumed 3.
# merge_blocks.R was written for this, but block_merge_map.csv is header-only.
#
# WHY IT WAS EMPTY: a join bug, not a data problem. It stripped the vote
# table's blk ("Drosera_regia_vs_Nepenthes_gracilis: 578") down to "578" and
# joined that against as.character(blkID) from the coordinate file. If those
# formats differ, every row drops and filter(!is.na(s)) empties the table.
# So merging was never actually tested.
#
# THIS SCRIPT tests it directly, geometry only. merge_blocks.R also required
# vote CONCORDANCE between neighbours, which is sensible for preserving
# exchange boundaries but adds a second way to fail. Geometry first: if that
# cannot reach ~3 units per region, concordance will not help either.
#
# MERGE RULE: same species, same ancestral region, same Drosera chromosome,
# adjacent in position, gap <= MAXGAP. Chain them.
#
# SUCCESS: median n_units per region drops toward 3 and top_frac rises.
# FAILURE: units stay >5 -> the vote test is dead at every available unit and
#          DR01 closes on the depth ratio alone.
#
# IN   DR/out/DR01_votes.csv, DR/locus_meta.tsv,
#      ../genespace/results/{syntenicBlock_coordinates.csv,combBed.txt}
# OUT  DR/out/DR01b_merge.csv
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2)})
setwd(Sys.getenv("SUBG_BASE", getwd()))
GSD <- file.path(dirname(getwd()),"genespace","results")
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

hr("1. block ID formats — the thing that broke merge_blocks.R")
blk <- read_csv(file.path(GSD,"syntenicBlock_coordinates.csv"), show_col_types=FALSE)
cat("  blkID class:", class(blk$blkID), "\n  examples:\n")
print(head(unique(blk$blkID), 6))
if (file.exists("tract_votes_blocks7.csv")) {
  v <- read_csv("tract_votes_blocks7.csv", show_col_types=FALSE, n_max=5000)
  cat("\n  vote-table blk examples:\n"); print(head(unique(v$blk), 3))
  cat("  after sub('.*:\\\\s*','',blk):\n")
  print(head(unique(sub(".*:\\s*", "", v$blk)), 6))
  cat("\n  do the stripped ids match blkID?  ",
      any(sub(".*:\\s*","",v$blk) %in% as.character(blk$blkID)), "\n")
}

hr("2. block extents on the Drosera side")
ext <- bind_rows(
  blk %>% filter(genome2=="Nepenthes_gracilis", genome1 %in% DROS) %>%
    transmute(sp=genome1, nep_chr=chr2, sp_chr=chr1, bid=as.character(blkID),
              s=pmin(startBp1,endBp1), e=pmax(startBp1,endBp1), nh=nHits1),
  blk %>% filter(genome1=="Nepenthes_gracilis", genome2 %in% DROS) %>%
    transmute(sp=genome2, nep_chr=chr1, sp_chr=chr2, bid=as.character(blkID),
              s=pmin(startBp2,endBp2), e=pmax(startBp2,endBp2), nh=nHits2)) %>%
  filter(!is.na(s), grepl("_dom$", nep_chr)) %>% distinct()
cat(sprintf("  blocks: %d across %d species\n", nrow(ext), n_distinct(ext$sp)))
cat("\n  gap to the previous block on the same chromosome (Mb):\n")
gaps <- ext %>% arrange(sp, sp_chr, s) %>% group_by(sp, sp_chr) %>%
  mutate(gap_mb=(s-lag(e))/1e6) %>% ungroup() %>% filter(!is.na(gap_mb))
print(summary(gaps$gap_mb))
cat(sprintf("  gaps <= 1 Mb: %.0f%% | <= 5 Mb: %.0f%%\n",
            100*mean(gaps$gap_mb<=1), 100*mean(gaps$gap_mb<=5)))
cat("  MAXGAP must be large enough to chain, small enough not to fuse\n")
cat("  everything into one unit per chromosome.\n")

hr("3. merge on geometry, sweeping MAXGAP")
merge_at <- function(G) {
  ext %>% arrange(sp, nep_chr, sp_chr, s) %>%
    group_by(sp, nep_chr, sp_chr) %>%
    mutate(newg = is.na(lag(e)) | (s - lag(e)) > G,
           m = cumsum(newg)) %>% ungroup() %>%
    mutate(munit = sprintf("%s_%s_m%03d", sp_chr, nep_chr, m))
}
for (G in c(0.5e6, 2e6, 5e6, 20e6)) {
  M <- merge_at(G)
  cat(sprintf("  MAXGAP %5.1f Mb : %5d blocks -> %5d merged units\n",
              G/1e6, nrow(M), n_distinct(paste(M$sp, M$munit))))
}

hr("4. does it fix n_units per region?")
votes <- read_csv("DR/out/DR01_votes.csv", show_col_types=FALSE)
LOC <- read_tsv("DR/locus_meta.tsv", show_col_types=FALSE) %>%
       group_by(tip) %>% slice(1) %>% ungroup()
bedpos <- read.table(file.path(GSD,"combBed.txt"), header=TRUE, sep="\t",
                     quote="", comment.char="", stringsAsFactors=FALSE) %>%
  filter(genome %in% DROS) %>% transmute(genome, gene=id, mid=(start+end)/2)
sel <- LOC %>% filter(genome %in% DROS) %>%
  group_by(locus, genome) %>% filter(n()==3) %>% ungroup() %>%
  mutate(region=sub("-.*$","",locus)) %>% left_join(bedpos, by=c("genome","gene"))

assign_unit <- function(M) {
  key <- M %>% select(sp, nep_chr, sp_chr, s, e, munit)
  f <- function(g, r, c, m) {
    if (is.na(m)) return(NA_character_)
    h <- key[key$sp==g & key$nep_chr==r & key$sp_chr==c & key$s<=m & key$e>=m, ]
    if (!nrow(h)) NA_character_ else h$munit[1]
  }
  mapply(f, sel$genome, sel$region, sel$chr, sel$mid)
}
res <- bind_rows(lapply(c(0.5e6, 2e6, 5e6, 20e6), function(G) {
  sel$unit <- assign_unit(merge_at(G))
  d <- sel %>% filter(!is.na(unit)) %>%
    group_by(locus, genome) %>% filter(n()==3) %>% ungroup()
  d %>% group_by(genome, region) %>%
    summarise(loci=n_distinct(locus), units=n_distinct(unit), .groups="drop") %>%
    mutate(maxgap_mb=G/1e6)
}))
cat("  units used per (species, region) at each MAXGAP:\n")
cat("  TARGET: 3. Anything above means the chi-square null of 1/3 is wrong.\n\n")
print(as.data.frame(res %>% group_by(genome, maxgap_mb) %>%
  summarise(regions=n(), median_units=median(units), min_units=min(units),
            max_units=max(units), .groups="drop") %>%
  pivot_wider(names_from=maxgap_mb, values_from=c(median_units,max_units),
              id_cols=c(genome, regions))), row.names=FALSE)
write_csv(res, "DR/out/DR01b_merge.csv")

hr("VERDICT")
best <- res %>% group_by(maxgap_mb) %>%
  summarise(median_units=median(units),
            frac_at_3=round(mean(units<=3),3),
            frac_le5=round(mean(units<=5),3), .groups="drop")
print(as.data.frame(best), row.names=FALSE)
cat("\n  If some MAXGAP gets median_units to ~3 with frac_at_3 well above 0,\n")
cat("  the vote test is rescuable and DR01 should be re-run with that unit.\n")
cat("  If units stay >5 everywhere, merging cannot fix it: the region simply\n")
cat("  contains more than three independent segments, and DR01 closes on the\n")
cat("  depth ratio alone.\n")

hr("5. REGIA VOTE TEST — merged units at MAXGAP = 2 Mb")
cat("  regia is the ONLY species reaching 3 merged units per region, so it is\n")
cat("  the only one where the 3-category chi-square null is valid. It is also\n")
cat("  flat from 2 Mb to 20 Mb, meaning its units are separated by >20 Mb --\n")
cat("  a real property, not an artefact of the threshold.\n")
cat("  Power at regia's ~64 loci/region: 0.98 at a true winning fraction of\n")
cat("  0.60, 0.70 at 0.50. Well powered.\n\n")
MG <- 2e6
Mr <- merge_at(MG)
selr <- sel %>% filter(genome=="Drosera_regia")
selr$unit <- local({
  key <- Mr %>% filter(sp=="Drosera_regia") %>% select(nep_chr, sp_chr, s, e, munit)
  mapply(function(r,c,m) {
    if (is.na(m)) return(NA_character_)
    h <- key[key$nep_chr==r & key$sp_chr==c & key$s<=m & key$e>=m, ]
    if (!nrow(h)) NA_character_ else h$munit[1]
  }, selr$region, selr$chr, selr$mid)
})
cat(sprintf("  regia copies assigned to a merged unit: %d of %d (%.1f%%)\n",
            sum(!is.na(selr$unit)), nrow(selr), 100*mean(!is.na(selr$unit))))
selr <- selr %>% filter(!is.na(unit)) %>%
        group_by(locus) %>% filter(n()==3, n_distinct(unit)==3) %>% ungroup()
cat(sprintf("  loci with all 3 copies in 3 DISTINCT merged units: %d\n",
            n_distinct(selr$locus)))
cat("  (loci where two copies land in one unit are excluded: the pair is\n")
cat("   ambiguous there, which is the same reason chromosomes failed.)\n")

VT <- votes %>% filter(genome=="Drosera_regia", clear) %>%
      select(locus, d_sis, d_out, ratio, stratum)
uu <- selr %>% select(locus, tip, unit)
tips3 <- uu %>% group_by(locus) %>% summarise(units=list(sort(unit)), .groups="drop")
K2 <- read_csv("DR/out/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS>=0, dS<5, codons>=100)
DD2 <- bind_rows(K2 %>% transmute(locus=anchor, a=seq1, b=seq2, d=dS),
                 K2 %>% transmute(locus=anchor, a=seq2, b=seq1, d=dS))
dmap2 <- setNames(DD2$d, paste(DD2$locus, DD2$a, DD2$b))
gd2 <- function(l,x,y) unname(dmap2[paste(l,x,y)])
nept2 <- LOC %>% filter(genome=="Nepenthes_gracilis") %>%
         group_by(locus) %>% slice(1) %>% ungroup() %>% transmute(locus, nep=tip)

RV <- bind_rows(lapply(unique(selr$locus), function(l) {
  d <- selr[selr$locus==l, ]
  N <- nept2$nep[nept2$locus==l][1]
  if (is.na(N) || nrow(d)!=3) return(NULL)
  t <- d$tip; u <- d$unit
  m <- matrix(NA_real_, 3, 3)
  for (x in 1:2) for (y in (x+1):3) { m[x,y] <- m[y,x] <- gd2(l,t[x],t[y]) }
  dn <- vapply(1:3, function(x) gd2(l,t[x],N), numeric(1))
  if (any(is.na(dn)) || any(is.na(m[upper.tri(m)]))) return(NULL)
  S <- c(m[1,2]+dn[3], m[1,3]+dn[2], m[2,3]+dn[1])
  pr <- list(c(1,2),c(1,3),c(2,3))[[which.min(S)]]
  og <- setdiff(1:3, pr); o <- sort(S)
  tibble(locus=l, region=sub("-.*$","",l),
         pair=paste(sort(u[pr]), collapse="|"),
         d_sis=m[pr[1],pr[2]], d_out=mean(c(m[pr[1],og], m[pr[2],og])),
         gap=(o[2]-o[1])/mean(S))
})) %>% filter(gap > 0.02)
cat(sprintf("\n  loci with a clear winner: %d\n", nrow(RV)))

chi3 <- function(v) { v <- c(v, rep(0L, max(0,3-length(v))))[1:3]
                      n <- sum(v); e <- n/3; sum((v-e)^2/e) }
REGV <- RV %>% group_by(region) %>%
  group_modify(function(x, key) {
    tb <- sort(table(x$pair), decreasing=TRUE)
    tibble(n=nrow(x), n_units=length(tb), top_pair=names(tb)[1],
           top_n=as.integer(tb)[1], chi2=chi3(as.integer(tb)))
  }) %>% ungroup() %>%
  mutate(top_frac=round(top_n/n,3), p=1-pchisq(chi2, df=2),
         p_adj=p.adjust(p,"BH"))
print(as.data.frame(REGV %>% filter(n>=15) %>%
  transmute(region, n, n_units, top_frac, chi2=round(chi2,2),
            p_adj=signif(p_adj,3))), row.names=FALSE)
cat(sprintf("\n  regions tested: %d | significant after BH: %d\n",
            sum(REGV$n>=15), sum(REGV$n>=15 & REGV$p_adj<0.05)))
cat(sprintf("  median top_frac: %.3f (uniform expectation 0.333)\n",
            median(REGV$top_frac[REGV$n>=15])))
cat("\n  CROSS-CHECK: is the winning pair's depth the ancestral A1/A2 split?\n")
win <- RV %>% inner_join(REGV %>% select(region, top_pair), by="region") %>%
       mutate(is_top = pair==top_pair)
print(as.data.frame(win %>% group_by(is_top) %>%
  summarise(loci=n(), med_d_sis=round(median(d_sis),3),
            med_d_out=round(median(d_out),3),
            ratio=round(median(d_sis)/median(d_out),3), .groups="drop")),
  row.names=FALSE)
cat("  the winning pair should sit near d_sis 0.22 (script 18 sec6). If the\n")
cat("  non-winning loci sit much lower, they are recent duplicates instead.\n")
write_csv(RV, "DR/out/DR01b_regia_votes.csv")
write_csv(REGV, "DR/out/DR01b_regia_regions.csv")

pv <- ggplot(RV %>% group_by(region) %>% mutate(nn=n()) %>% filter(nn>=15) %>% ungroup(),
             aes(region, fill=pair)) +
  geom_bar(position="fill") +
  geom_hline(yintercept=c(1/3,2/3), colour="grey30", linetype="dotted") +
  labs(title="DR01b - regia: which merged-unit pair is sister, per region",
       subtitle="merged at 2 Mb, 3 units per region | AAB = one colour dominates | dotted = 1/3, 2/3",
       x="ancestral region", y="fraction of loci", fill="unit pair") +
  theme_minimal(10) + theme(legend.position="bottom",
                            legend.text=element_text(size=6))
suppressWarnings(ggsave("DR/fig/DR01b_regia_votes.pdf", pv, width=10, height=7))
cat("\n  FIG: DR/fig/DR01b_regia_votes.pdf\n")
