#!/usr/bin/env Rscript
# sCF has comment lines before the header; read.table needs comment.char="#".
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
setwd(Sys.getenv("SUBG_BASE", getwd()))
for (M in c("A_full","B_full","all13_full")) {
  f <- sprintf("DR/tree/concat/%s_scf.cf.stat", M)
  if (!file.exists(f)) next
  tab <- read.table(f, header=TRUE, comment.char="#")
  cat(sprintf("\n=== %s ===\n", M))
  print(tab[, intersect(c("ID","sCF","sDF1","sDF2","sN","Label"), names(tab))])
  cat("  sCF is the %% of SITES supporting the branch. 33%% = no signal,\n")
  cat("  the three resolutions being equally supported.\n")
}
