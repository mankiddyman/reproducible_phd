#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); NBOOT <- 200
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
NEP <- "Nepenthes_gracilis"

emfit <- function(x, k, iters=800) {
  n <- length(x); mu <- as.numeric(quantile(x, seq(.5/k, 1-.5/k, length.out=k)))
  sdv <- rep(sd(x)/k, k); w <- rep(1/k, k); ll_old <- -Inf; ll <- -Inf
  for (it in seq_len(iters)) {
    dens <- sapply(seq_len(k), function(j) w[j]*dnorm(x, mu[j], sdv[j]))
    if (k==1) dens <- matrix(dens, ncol=1)
    tot <- rowSums(dens); tot[tot < 1e-300] <- 1e-300; r <- dens/tot
    w <- colMeans(r); cs <- colSums(r)
    mu <- colSums(r*x)/cs
    sdv <- sqrt(colSums(r*(outer(x, mu, "-")^2))/cs); sdv[sdv < 1e-3] <- 1e-3
    ll <- sum(log(tot)); if (abs(ll-ll_old) < 1e-9) break; ll_old <- ll
  }
  o <- order(mu)
  list(mu=mu[o], sd=sdv[o], w=w[o], ll=ll, k=k, bic=-2*ll + (3*k-1)*log(n))
}

k <- read_csv("ks/pairwise_ks.csv", show_col_types=FALSE) %>%
     filter(!is.na(dS), dS>0, dS<5, codons>=100)
tm <- read_tsv("wgd7/tip_meta.tsv", show_col_types=FALSE)
key <- tm %>% select(anchor=nep_gene, k1=tip, genome, chr) %>% distinct()
WI <- k %>% filter(sp1==sp2, sp1!=NEP) %>% transmute(anchor, sp=sp1, g1=seq1, g2=seq2, dS) %>%
  left_join(key %>% rename(g1=k1, c1=chr), by=c("anchor","g1","sp"="genome")) %>%
  left_join(key %>% rename(g2=k1, c2=chr), by=c("anchor","g2","sp"="genome")) %>%
  filter(!is.na(c1), !is.na(c2), c1 != c2)

hr("1. mixture on log(dS), different-chromosome within-species pairs")
cat("  AAB -> 2 components, weights ~1:2.   ABC -> 1 component.\n\n")
set.seed(1)
fits <- lapply(split(WI, WI$sp), function(d) {
  x <- log(d$dS); best <- NULL
  for (kk in 1:4) { f <- emfit(x, kk); if (is.null(best) || f$bic < best$bic) best <- f }
  list(sp=d$sp[1], n=length(x), best=best,
       bics=vapply(1:4, function(kk) emfit(x,kk)$bic, numeric(1)))
})
print(as.data.frame(bind_rows(lapply(fits, function(f) tibble(
  sp=f$sp, n=f$n, k=f$best$k,
  bic1=round(f$bics[1]), bic2=round(f$bics[2]),
  bic3=round(f$bics[3]), bic4=round(f$bics[4]))))), row.names=FALSE)

hr("2. components (back-transformed to dS)")
comp <- bind_rows(lapply(fits, function(f) tibble(
  sp=f$sp, comp=seq_len(f$best$k), dS=round(exp(f$best$mu),3),
  weight=round(f$best$w,3), log_sd=round(f$best$sd,3))))
print(as.data.frame(comp), row.names=FALSE)

hr("3. CALIBRATION — does it recover regia?")
cat("  regia is AAB from four independent lines: A1|A2 = 0.22, A|B = 0.355,\n")
cat("  ratio 0.62, expected weights 1:2.\n\n")
rg <- comp %>% filter(sp=="Drosera_regia")
print(as.data.frame(rg), row.names=FALSE)
if (nrow(rg) >= 2) {
  cat(sprintf("\n  observed shallow/deep ratio: %.3f  (expected ~0.62)\n",
              rg$dS[1]/rg$dS[2]))
  cat(sprintf("  observed weight of shallow component: %.3f  (expected ~0.33)\n", rg$weight[1]))
  cat("  Close on both -> the method resolves AAB and the other rows are readable.\n")
} else cat("\n  Only 1 component for regia -> the method CANNOT resolve a known AAB.\n")

hr("4. bootstrap CIs on the chosen model")
bs <- bind_rows(lapply(split(WI, WI$sp), function(d) {
  x <- log(d$dS); kk <- fits[[d$sp[1]]]$best$k; if (kk < 2) return(NULL)
  m <- vapply(seq_len(NBOOT), function(i) {
    f <- tryCatch(emfit(sample(x, length(x), replace=TRUE), kk), error=function(e) NULL)
    if (is.null(f)) rep(NA_real_, 2) else c(exp(f$mu[1]), f$w[1])
  }, numeric(2))
  tibble(sp=d$sp[1], k=kk,
         c1_lo=round(quantile(m[1,], .025, na.rm=TRUE),3),
         c1_hi=round(quantile(m[1,], .975, na.rm=TRUE),3),
         w1_lo=round(quantile(m[2,], .025, na.rm=TRUE),3),
         w1_hi=round(quantile(m[2,], .975, na.rm=TRUE),3))
}))
if (nrow(bs)) print(as.data.frame(bs), row.names=FALSE)
cat("\n  Weight CI excluding 0.33 argues against a clean 1:2 AAB.\n")

hr("5. validation against the 3-copy loci")
cat("  At 3-chromosome loci the pairing is known. Do the shallow-component\n")
cat("  pairs correspond to the four-point sister pairs?\n\n")
tri <- read_csv("trees/qc/triplet_all_species.csv", show_col_types=FALSE)
if (nrow(tri)) print(as.data.frame(tri %>% group_by(genome) %>%
  summarise(triplets=n(), d_sis=round(median(d_sis),3),
            d_out=round(median(d_out),3), .groups="drop")), row.names=FALSE)
cat("\n  d_sis should track the shallow component, d_out the deep one.\n")

p <- WI %>% mutate(sp=sub("Drosera_","D_",sp)) %>%
  ggplot(aes(dS)) + geom_histogram(bins=60, fill="#378ADD", colour="white") +
  facet_wrap(~sp, scales="free_y") + scale_x_log10() +
  labs(title="FIG52 - within-species dS, different chromosomes",
       subtitle="two peaks at a 1:2 weight ratio = AAB; one peak = ABC",
       x="dS (log scale)", y="pairs") + theme_minimal(10)
suppressWarnings(ggsave("trees/qc/FIG52_mixture.pdf", p, width=10, height=6))
write_csv(comp, "trees/qc/mixture_components.csv")
