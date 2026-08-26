#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Each Dionaea pair gives a LOCAL sgA/sgB with an arbitrary label. A Drosera block that
# votes on TWO different pairs constrains those two labels: same side => same subgenome.
# Build the signed graph, two-colour it, count violated constraints against a permuted null.
# Uses NO fractionation information.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
set.seed(1); MINV <- 5; MINF <- 0.70; B <- 999

v <- read_csv("tract_votes_blocks7.csv", show_col_types = FALSE) %>%
  mutate(species = sub("Drosera_","",sp), bid = sub(".*:\\s*","",blk),
         dblk = sprintf("%s_%s_b%s", substr(species,1,3),
                        sub("_hap1$|_collapsed$","",lin_chr), bid))

## A GENESPACE block is defined against ONE Nepenthes chromosome, which maps to ONE
## Dionaea pair -- so a block can never link two pairs. The linking unit must be the
## Drosera CHROMOSOME, which carries blocks against several Nepenthes chromosomes.
v <- v %>% mutate(dchr = sprintf("%s_%s", substr(species,1,3),
                                 sub("_hap1$|_collapsed$","",lin_chr)))
C <- v %>% group_by(dchr, species, pair) %>%
  summarise(n = n(), fA = mean(vote == "A"), .groups = "drop") %>%
  filter(n >= MINV, fA >= MINF | fA <= 1 - MINF) %>%
  mutate(side = ifelse(fA >= MINF, 1L, -1L), dblk = dchr)
cat(sprintf("confident (chromosome, pair) calls: %d | Drosera chromosomes: %d\n",
            nrow(C), n_distinct(C$dchr)))
cat("\ncalls per Dionaea pair:\n"); print(as.data.frame(count(C, pair)))

link <- C %>% group_by(dchr) %>% filter(n_distinct(pair) >= 2) %>% ungroup()
cat(sprintf("\nDrosera chromosomes calling on >=2 Dionaea pairs: %d\n", n_distinct(link$dchr)))
if (n_distinct(link$dchr)) {
  cat("\nlinking chromosomes and the pairs they span:\n")
  print(as.data.frame(link %>% group_by(dchr) %>%
    summarise(pairs = paste(sort(pair), collapse=","), sides = paste(side, collapse=","),
              votes = sum(n), .groups="drop") %>% arrange(desc(votes))), right = FALSE)
}
if (!n_distinct(link$dchr)) { cat("still no cross-pair links\n"); quit(save="no") }

edges <- function(L) bind_rows(lapply(split(L, L$dchr), function(d) {
  if (nrow(d) < 2) return(NULL)
  cb <- combn(seq_len(nrow(d)), 2)
  tibble(dblk = d$dchr[1], species = d$species[1],
         p1 = d$pair[cb[1,]], p2 = d$pair[cb[2,]],
         same = d$side[cb[1,]] == d$side[cb[2,]]) }))
E <- edges(link)
cat(sprintf("constraints: %d (%d SAME, %d DIFFERENT)\n", nrow(E), sum(E$same), sum(!E$same)))
cat("\nconstraints per species:\n"); print(as.data.frame(count(E, species)))
cat("\nconstraint matrix (pair x pair, SAME/total):\n")
print(as.data.frame(E %>% mutate(a=pmin(p1,p2), b=pmax(p1,p2)) %>%
  group_by(a,b) %>% summarise(n=n(), same=sum(same), .groups="drop") %>%
  arrange(desc(n))))

colour2 <- function(E) {
  nodes <- sort(unique(c(E$p1, E$p2)))
  col <- setNames(rep(NA_integer_, length(nodes)), nodes)
  nb <- lapply(nodes, function(x) which(E$p1 == x | E$p2 == x)); names(nb) <- nodes
  for (s in nodes) {
    if (!is.na(col[s])) next
    col[s] <- 1L; q <- s
    while (length(q)) { x <- q[1]; q <- q[-1]
      for (i in nb[[x]]) { y <- if (E$p1[i] == x) E$p2[i] else E$p1[i]
        w <- if (E$same[i]) col[x] else 3L - col[x]
        if (is.na(col[y])) { col[y] <- w; q <- c(q, y) } } }
  }
  viol <- sum(mapply(function(a,b,s) (col[a] == col[b]) != s, E$p1, E$p2, E$same))
  list(col = col, viol = viol)
}
r <- colour2(E)
nul <- replicate(B, colour2(mutate(E, same = sample(same)))$viol)
p <- (1 + sum(nul <= r$viol)) / (B + 1)
cat(sprintf("\n=== two-colouring: %d/%d constraints violated (null median %.0f, p = %.4f) ===\n",
            r$viol, nrow(E), median(nul), p))
print(data.frame(pair = names(r$col), group = ifelse(r$col == 1, "A", "B")))

## agreement with the fractionation partition -- reported, never used above
ref <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE) %>%
  filter(anchor == "Nepenthes_gracilis") %>%
  transmute(pair = pair_lab, X = winner, chrA = pmin(chrA, chrB))
cmp <- v %>% distinct(pair, chrA) %>% left_join(ref, by = "pair") %>%
  mutate(sgA_is_X = chrA.x == X,
         topo = ifelse(r$col[pair] == 1, "A", "B")) %>%
  filter(!is.na(sgA_is_X))
cat("\n=== independent check: does the topological grouping match fractionation? ===\n")
print(as.data.frame(cmp %>% select(pair, sgA_is_X, topo)))
k <- sum((cmp$topo == "A") == cmp$sgA_is_X); n <- nrow(cmp)
cat(sprintf("agree %d/%d | sign test p = %.4f (floor with one label fixed = %.4f)\n",
            max(k, n-k), n, binom.test(max(k, n-k), n, 0.5)$p.value, 0.5^(n-1)))
write_csv(E, "global_constraints.csv")
