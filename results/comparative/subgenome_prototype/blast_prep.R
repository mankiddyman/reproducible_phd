#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Query sets for the sequence-level presence test.
#   lost : survivor-side protein, searched against the homeolog that "lost" it
#   pos  : protein from a 1:2 retained HOG -> homolog IS present and annotated
#   neg  : same queries against a NON-homeologous chromosome (background)
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
set.seed(1); N <- 200

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
res <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE)

part <- res %>% filter(anchor == "Nepenthes_gracilis") %>%
  transmute(pair_lab, X = winner, Y = ifelse(winner == chrA, chrB, chrA))

b <- bed %>% mutate(isRep = as.logical(isArrayRep)) %>% filter(is.na(isRep) | isRep)
anchored <- b %>% filter(genome != "Dionaea_muscipula") %>% distinct(globHOG) %>% pull(globHOG)

dio <- b %>% filter(genome == "Dionaea_muscipula") %>%
  mutate(pl = sub("_sg[0-9]+_s[0-9]+$", "", chr)) %>%
  filter(globHOG %in% anchored) %>%
  select(gene = id, chr, pl, globHOG)

st <- dio %>% group_by(globHOG) %>%
  summarise(n_dio = n(), n_pl = n_distinct(pl), n_chr = n_distinct(chr), .groups = "drop")
dio <- left_join(dio, st, by = "globHOG")

samp <- function(d) if (nrow(d) > N) d[sample(nrow(d), N), ] else d

jobs <- list()
for (i in seq_len(nrow(part))) {
  p <- part$pair_lab[i]; X <- part$X[i]; Y <- part$Y[i]
  negT <- part$Y[if (i == nrow(part)) 1 else i + 1]   # chromosome from another pair

  lostY <- samp(filter(dio, n_dio == 1, chr == X))          # X kept it, Y "lost" it
  lostX <- samp(filter(dio, n_dio == 1, chr == Y))
  posX  <- samp(filter(dio, n_dio == 2, n_pl == 1, n_chr == 2, chr == X, pl == p))
  posY  <- samp(filter(dio, n_dio == 2, n_pl == 1, n_chr == 2, chr == Y, pl == p))

  for (j in list(list("lostY", lostY, Y), list("lostX", lostX, X),
                 list("posX",  posX,  Y), list("posY",  posY,  X),
                 list("negX",  lostY, negT))) {
    if (!nrow(j[[2]])) next
    tag <- sprintf("%s_%s", p, j[[1]])
    writeLines(j[[2]]$gene, file.path("blast/prot", paste0(tag, ".ids")))
    jobs[[length(jobs)+1]] <- tibble(tag = tag, pair_lab = p, cat = j[[1]],
                                     target = j[[3]], n = nrow(j[[2]]))
  }
}
jb <- bind_rows(jobs)
write_csv(jb, "blast/jobs.csv")
write_csv(part, "blast/partition.csv")
writeLines(unique(c(part$X, part$Y)), "blast/chroms.txt")
cat("=== query sets ===\n"); print(as.data.frame(jb))
