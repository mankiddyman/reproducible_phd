#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(ape); library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); NBOOT <- 300
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"; DIO <- "Dionaea_muscipula"
DROS <- c("Drosera_regia","Drosera_binata","Drosera_paradoxa",
          "Drosera_scorpioides","Drosera_capensis")
ALL <- c(DROS, DIO, NEP)

k  <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
      filter(!is.na(dS), dS>0, dS<3, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
sc <- tm %>% count(nep_gene, genome, name="k") %>% filter(k==1)

## Drosera + Nepenthes: true single-copy orthologs
SP <- k %>% filter(sp1 != sp2) %>%
  semi_join(sc, by=c("anchor"="nep_gene","sp1"="genome")) %>%
  semi_join(sc, by=c("anchor"="nep_gene","sp2"="genome")) %>%
  mutate(a=pmin(sp1,sp2), b=pmax(sp1,sp2)) %>% select(anchor,a,b,dS)

## Dionaea is 2-copy by the wgd7 filter, so take the nearest copy per locus
## (the ortholog under shared A/B) and also the mean, to bracket the bias.
DI <- k %>% filter(xor(sp1==DIO, sp2==DIO)) %>%
  mutate(o=ifelse(sp1==DIO, sp2, sp1)) %>% filter(o != DIO) %>%
  semi_join(sc, by=c("anchor"="nep_gene","o"="genome")) %>%
  group_by(anchor, o) %>%
  summarise(dmin=min(dS), dmean=mean(dS), n=n(), .groups="drop") %>% filter(n==2)

hr("0. distance matrix")
M <- SP %>% group_by(a,b) %>% summarise(d=median(dS), n=n(), .groups="drop")
Dm <- DI %>% group_by(o) %>%
  summarise(n=n(), dmin=round(median(dmin),3), dmean=round(median(dmean),3), .groups="drop")
cat("  Drosera/Nepenthes single-copy pairs:\n")
print(as.data.frame(M %>% mutate(across(c(a,b), ~sub("Drosera_","",sub("Nepenthes_gracilis","NEP",.))))), row.names=FALSE)
cat("\n  Dionaea (2-copy): nearest vs mean copy\n")
print(as.data.frame(Dm %>% mutate(o=sub("Drosera_","",sub("Nepenthes_gracilis","NEP",o)))), row.names=FALSE)
cat("\n  min is biased low (min of 2 draws), mean biased high. Both fitted below.\n")

fit_rates <- function(Mx, Dx, which=c("dmin","dmean")) {
  which <- match.arg(which)
  mat <- matrix(NA_real_, length(ALL), length(ALL), dimnames=list(ALL, ALL))
  diag(mat) <- 0
  for (i in seq_len(nrow(Mx))) { mat[Mx$a[i], Mx$b[i]] <- Mx$d[i]; mat[Mx$b[i], Mx$a[i]] <- Mx$d[i] }
  for (i in seq_len(nrow(Dx))) { mat[DIO, Dx$o[i]] <- Dx[[which]][i]; mat[Dx$o[i], DIO] <- Dx[[which]][i] }
  if (any(is.na(mat))) return(NULL)
  tr <- nj(as.dist(mat)); tr <- root(tr, outgroup=NEP, resolve.root=TRUE)
  r2t <- node.depth.edgelength(tr)[seq_len(Ntip(tr))]; names(r2t) <- tr$tip.label
  r2t[ALL]
}
Dx0 <- DI %>% group_by(o) %>% summarise(dmin=median(dmin), dmean=median(dmean), .groups="drop")

hr("1. root-to-tip from a fitted NJ tree = relative rate")
cat("  Tips are contemporaneous, so root-to-tip is rate x total time.\n\n")
out <- bind_rows(lapply(c("dmin","dmean"), function(w) {
  r <- fit_rates(M, Dx0, w); if (is.null(r)) return(NULL)
  tibble(dio_est=w, species=names(r), r2t=round(r,4),
         rate=round(r/r[["Drosera_binata"]],3))
}))
print(as.data.frame(out %>% select(species, dio_est, r2t, rate) %>%
  pivot_wider(names_from=dio_est, values_from=c(r2t, rate)) %>%
  arrange(desc(rate_dmin))), row.names=FALSE)

set.seed(1)
bs <- vapply(seq_len(NBOOT), function(i) {
  Mb <- SP %>% group_by(a,b) %>% summarise(d=median(sample(dS,n(),replace=TRUE)), .groups="drop")
  Db <- DI %>% group_by(o) %>% summarise(dmin=median(sample(dmin,n(),replace=TRUE)),
                                         dmean=median(dmin), .groups="drop")
  r <- tryCatch(fit_rates(Mb, Db, "dmin"), error=function(e) NULL)
  if (is.null(r)) rep(NA_real_, length(ALL)) else r/r[["Drosera_binata"]]
}, numeric(length(ALL)))
cat("\n  bootstrap 95% CI on rate (dmin fit):\n")
print(as.data.frame(tibble(species=ALL,
  lo=round(apply(bs,1,quantile,.025,na.rm=TRUE),3),
  hi=round(apply(bs,1,quantile,.975,na.rm=TRUE),3))), row.names=FALSE)
cat("\n  Ad-hoc model predicting d1(Dio,regia)=0.278 vs 0.276 observed implied\n")
cat("  regia 0.256, Dionaea 0.400 relative to binata. Compare.\n")

RT <- out %>% filter(dio_est=="dmin") %>% select(sp=species, rate)

hr("2. cross-check against distance to Nepenthes")
print(as.data.frame(bind_rows(
  M %>% filter(a==NEP|b==NEP) %>% transmute(sp=ifelse(a==NEP,b,a), to_nep=round(d,3)),
  Dx0 %>% filter(o==NEP) %>% transmute(sp=DIO, to_nep=round(dmin,3))) %>%
  left_join(RT, by="sp") %>% arrange(desc(rate))), row.names=FALSE)
cat("\n  FIG41 reported ~0.85 (regia) to ~1.4 (scorpioides).\n")

hr("3. THE TEST — is the A/B split one shared event?")
cat("  T = dS(cross-AB pairs) / (2 x rate). If A/B is one split, these must\n")
cat("  be EQUAL across species. Raw values span 0.358 to 1.416.\n\n")
key <- tm %>% select(anchor=nep_gene, k1=tip, genome, chr) %>% distinct() %>%
       mutate(region=sub("-.*$","",anchor))
seg <- read_csv("trees/qc/segment_delta_calibrated.csv", show_col_types=FALSE) %>%
  mutate(side=case_when(med_delta < -0.08 ~ "A", med_delta > 0.08 ~ "B", TRUE ~ "amb")) %>%
  select(genome, region, chr, side)
P <- k %>% filter(sp1==sp2, sp1!=NEP) %>% transmute(anchor, sp=sp1, g1=seq1, g2=seq2, dS) %>%
  left_join(key %>% select(anchor, g1=k1, c1=chr, region), by=c("anchor","g1")) %>%
  left_join(key %>% select(anchor, g2=k1, c2=chr), by=c("anchor","g2")) %>%
  filter(!is.na(c1), !is.na(c2), c1!=c2) %>%
  left_join(seg %>% rename(sp=genome, c1=chr, s1=side), by=c("sp","region","c1")) %>%
  left_join(seg %>% rename(sp=genome, c2=chr, s2=side), by=c("sp","region","c2")) %>%
  filter(!is.na(s1), !is.na(s2), s1!="amb", s2!="amb") %>%
  mutate(cls=ifelse(s1==s2, paste0("same-",s1), "cross-AB")) %>%
  left_join(RT, by="sp")
print(as.data.frame(P %>% group_by(sp, cls) %>%
  summarise(pairs=n(), raw=round(median(dS),3), rate=round(rate[1],3),
            T=round(median(dS)/(2*rate[1]),3), .groups="drop") %>%
  arrange(cls, T)), row.names=FALSE)
cat("\n  cross-AB T collapsing onto one value = shared A/B, strongest form.\n")
cat("  same-A T = when each lineage's second A subgenome arose.\n")

hr("4. rate-corrected speciation depths — one axis with section 3")
print(as.data.frame(M %>% filter(a!=NEP, b!=NEP) %>%
  left_join(RT %>% rename(a=sp, ra=rate), by="a") %>%
  left_join(RT %>% rename(b=sp, rb=rate), by="b") %>%
  mutate(T=round(d/(ra+rb),3),
         pair=paste(sub("Drosera_","",a), sub("Drosera_","",b), sep="/")) %>%
  select(pair, dS=d, T) %>% arrange(T)), row.names=FALSE)
write_csv(out, "trees/qc/lineage_rates.csv")
