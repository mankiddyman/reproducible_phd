#!/usr/bin/env Rscript
# ESS + convergence for every MCMCtree run under DR/dating/.
# Read-only. Prints a table; writes DR/out/DRess.csv
suppressMessages(library(data.table)); options(width=200)
setwd(Sys.getenv("SUBG_BASE", getwd()))
ess1 <- function(x){ x <- x[is.finite(x)]; n <- length(x)
  if (n < 20 || var(x) == 0) return(NA_real_)
  a <- acf(x, lag.max=min(2000, n-1), plot=FALSE)$acf[-1]
  k <- which(a < 0.05)[1]; if (is.na(k)) k <- length(a)
  n/(1 + 2*sum(a[1:k])) }
rows <- list()
for (d in list.dirs("DR/dating", recursive=FALSE)) {
  f <- file.path(d, "mcmc.txt"); if (!file.exists(f)) next
  m <- tryCatch(fread(f), error=function(e) NULL); if (is.null(m) || !nrow(m)) next
  tc <- grep("^t_n", names(m), value=TRUE); if (!length(tc)) tc <- setdiff(names(m), "Gen")[1]
  e <- vapply(tc, function(cn) ess1(m[[cn]]), numeric(1))
  h <- m[[tc[1]]]; n <- length(h)
  rows[[length(rows)+1]] <- data.table(
    run = basename(d), samples = n, n_nodes = length(tc),
    min_ESS = round(min(e, na.rm=TRUE)), median_ESS = round(median(e, na.rm=TRUE)),
    worst_node = tc[which.min(e)],
    drift_pct = round(100*(mean(tail(h, n%/%2)) - mean(head(h, n%/%2)))/mean(h), 2),
    lnL_ESS = if ("lnL" %in% names(m)) round(ess1(m$lnL)) else NA_real_)
}
R <- rbindlist(rows)[order(min_ESS)]
R[, verdict := fifelse(min_ESS >= 200, "OK", fifelse(min_ESS >= 100, "MARGINAL", "FAIL"))]
print(R, row.names=FALSE)
fwrite(R, "DR/out/DRess.csv")
cat("\n  min ESS >= 200 is the usual bar. drift_pct compares the mean of the\n")
cat("  second half against the first: a few percent is fine, more is not mixed.\n")
cat("\n  runs failing the bar:", paste(R$run[R$min_ESS < 200], collapse=", "), "\n")
