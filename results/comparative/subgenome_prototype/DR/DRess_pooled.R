#!/usr/bin/env Rscript
# Pooled ESS + between-chain agreement for paired MCMCtree runs.
suppressMessages(library(data.table)); options(width=200)
setwd(Sys.getenv("SUBG_BASE", getwd()))
ess1 <- function(x){ x <- x[is.finite(x)]; n <- length(x)
  if (n < 20 || var(x) == 0) return(NA_real_)
  a <- acf(x, lag.max=min(2000, n-1), plot=FALSE)$acf[-1]
  k <- which(a < 0.05)[1]; if (is.na(k)) k <- length(a)
  n/(1 + 2*sum(a[1:k])) }
out <- list()
for (p in c("all13_ft","A_ft","B_ft","A_long","B_long")) {
  f1 <- sprintf("DR/dating/%s1/mcmc.txt", p); f2 <- sprintf("DR/dating/%s2/mcmc.txt", p)
  if (!file.exists(f1) || !file.exists(f2)) next
  m1 <- fread(f1); m2 <- fread(f2)
  tc <- intersect(grep("^t_n", names(m1), value=TRUE), grep("^t_n", names(m2), value=TRUE))
  e1 <- vapply(tc, function(c) ess1(m1[[c]]), numeric(1))
  e2 <- vapply(tc, function(c) ess1(m2[[c]]), numeric(1))
  mu1 <- vapply(tc, function(c) mean(m1[[c]]), numeric(1))
  mu2 <- vapply(tc, function(c) mean(m2[[c]]), numeric(1))
  sd1 <- vapply(tc, function(c) sd(m1[[c]]),  numeric(1))
  out[[p]] <- data.table(pair=p, nodes=length(tc),
    min_ESS_1=round(min(e1)), min_ESS_2=round(min(e2)),
    min_ESS_pooled=round(min(e1+e2)),
    max_mean_diff_pct=round(100*max(abs(mu1-mu2)/mu1),2),
    max_diff_in_sd=round(max(abs(mu1-mu2)/sd1),3),
    verdict=fifelse(min(e1+e2) >= 200, "OK pooled", "STILL SHORT"))
}
R <- rbindlist(out); print(R, row.names=FALSE)
cat("\n  max_diff_in_sd < 0.1 means the two chains agree to well within their own\n")
cat("  posterior width -- the standard MCMCtree convergence check. If that holds,\n")
cat("  pooling is legitimate and min_ESS_pooled is the number to report.\n")
