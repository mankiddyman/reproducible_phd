#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# Independently recompute every number the figure displays, from tract_votes_blocks7m.csv,
# and diff against window_matrix_rows.csv. Any mismatch is a bug in the plotting path.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
W <- 8; MINROW <- 6
ABB <- c(regia="reg", binata="bin", paradoxa="par", scorpioides="sco", capensis="cap")
FAIL <- 0L
chk <- function(ok, msg, detail = NULL) {
  cat(sprintf("[%s] %s\n", ifelse(ok, " OK ", "FAIL"), msg))
  if (!ok) { FAIL <<- FAIL + 1L; if (!is.null(detail)) print(as.data.frame(head(detail, 12))) }
}

v <- read_csv("tract_votes_blocks7m.csv", show_col_types = FALSE) %>%
  mutate(species = sub("Drosera_","",sp), bid = sub(".*:\\s*","",blk),
         unit = sprintf("%s %s b%s", ABB[species],
                        sub("_hap1$|_collapsed$","", sub("^chr","c",lin_chr)), bid))
R <- read_csv("window_matrix_rows.csv", show_col_types = FALSE)
wi <- read_csv("window_index2.csv", show_col_types = FALSE)

cat("=== A. inputs ===\n")
chk(all(v$vote %in% c("A","B")), "vote column is only A/B")
chk(!any(is.na(v$posA)), "no missing Dionaea positions")
chk(nrow(distinct(v, pair, anchor, sp_gene)) == nrow(v),
    "one row per (pair, anchor, voting gene) - no duplicated votes",
    v %>% count(pair, anchor, sp_gene) %>% filter(n > 1))
d2 <- v %>% distinct(pair, chrA, chrB) %>% count(pair) %>% filter(n > 1)
chk(nrow(d2) == 0, "each Dionaea pair maps to exactly one (chrA, chrB)", d2)

cat("\n=== B. windowing ===\n")
an <- v %>% distinct(pair, anchor, posA) %>% arrange(pair, posA) %>%
  group_by(pair) %>% mutate(win2 = (row_number()-1) %/% W + 1) %>% ungroup()
chk(nrow(an) == n_distinct(paste(v$pair, v$anchor)),
    "anchors are unique within a pair (a duplicate would shift all windows)")
w2 <- an %>% group_by(pair, win2) %>%
  summarise(n_anch2 = n(), mid2 = median(posA)/1e6, .groups="drop")
cmpw <- wi %>% select(pair, win, n_anch, mid_Mb) %>%
  full_join(w2, by = c("pair","win"="win2")) %>%
  mutate(bad = is.na(n_anch) | is.na(n_anch2) | n_anch != n_anch2 |
               abs(mid_Mb - mid2) > 1e-6)
chk(!any(cmpw$bad), "window index matches an independent recompute",
    filter(cmpw, bad))
mono <- an %>% group_by(pair) %>% summarise(ok = all(diff(posA) >= 0), .groups="drop")
chk(all(mono$ok), "anchors within a pair are position-sorted", filter(mono, !ok))
sz <- wi %>% group_by(pair) %>% summarise(bad = any(n_anch > W), .groups="drop")
chk(!any(sz$bad), sprintf("no window holds more than W=%d anchors", W), filter(sz, bad))

cat("\n=== C. block rows ===\n")
vw <- v %>% left_join(select(an, pair, anchor, win = win2), by = c("pair","anchor"))
B2 <- vw %>% group_by(pair, win, species, unit) %>%
  summarise(nv2 = n(), fA2 = mean(vote == "A"), .groups="drop")
keep2 <- B2 %>% group_by(pair, unit) %>%
  filter(sum(nv2) >= MINROW, n_distinct(win) >= 2) %>% ungroup() %>% distinct(pair, unit)
B2k <- semi_join(B2, keep2, by = c("pair","unit"))
Rb <- filter(R, kind == "block")
cmpb <- full_join(Rb %>% select(pair, win, species, unit, nv, fA),
                  B2k, by = c("pair","win","species","unit")) %>%
  mutate(bad = is.na(nv) | is.na(nv2) | nv != nv2 | abs(fA - fA2) > 1e-9)
chk(!any(cmpb$bad), sprintf("all %d block cells match an independent recompute", nrow(B2k)),
    filter(cmpb, bad) %>% select(pair, win, unit, nv, nv2, fA, fA2))
chk(setequal(paste(Rb$pair, Rb$unit), paste(keep2$pair, keep2$unit)),
    "the displayed block set equals the MINROW filter applied independently")

cat("\n=== D. species consensus ===\n")
S2 <- B2k %>% group_by(pair, win, species) %>%
  summarise(fA2 = sum(fA2 * nv2) / sum(nv2), nv2 = sum(nv2), .groups="drop")
Rs <- filter(R, kind == "spcons") %>% mutate(species2 = sub("^~ (\\w+) consensus$","\\1", unit))
Rs$species2 <- names(ABB)[match(Rs$species2, ABB)]
cmps <- full_join(Rs %>% select(pair, win, species, nv, fA), S2,
                  by = c("pair","win","species")) %>%
  mutate(bad = is.na(nv) | is.na(nv2) | nv != nv2 | abs(fA - fA2) > 1e-9)
chk(!any(cmps$bad), "each species consensus == pooled displayed blocks of that species",
    filter(cmps, bad) %>% select(pair, win, species, nv, nv2, fA, fA2))

cat("\n=== E. ALL consensus ===\n")
A2 <- S2 %>% group_by(pair, win) %>%
  summarise(nsp2 = n(), fA2 = mean(fA2), .groups = "drop")
Ra <- filter(R, kind == "allcons")
cmpa <- full_join(Ra %>% select(pair, win, nv, fA), A2, by = c("pair","win")) %>%
  mutate(bad = is.na(fA) | is.na(fA2) | abs(fA - fA2) > 1e-9 | nv != nsp2)
chk(!any(cmpa$bad), "ALL row == unweighted mean of species rows, nv == number of species",
    filter(cmpa, bad) %>% select(pair, win, nv, nsp2, fA, fA2))
chk(all(Ra$nv <= 5), "ALL row nv never exceeds 5 species", filter(Ra, nv > 5))

cat("\n=== F. ranges and totals ===\n")
chk(all(R$fA >= 0 & R$fA <= 1), "all fractions in [0,1]", filter(R, fA < 0 | fA > 1))
chk(all(R$nv > 0), "all counts positive", filter(R, nv <= 0))
tot_disp <- sum(Rb$nv)
chk(tot_disp == sum(B2k$nv2), sprintf("displayed vote total %d matches recompute", tot_disp))
cat(sprintf("     raw votes %d | displayed %d (%.0f%%) | dropped by MINROW %d\n",
            nrow(v), tot_disp, 100*tot_disp/nrow(v), nrow(v) - tot_disp))

cat("\n=== G. per-page summary ===\n")
pg <- Rb %>% group_by(pair) %>%
  summarise(rows = n_distinct(unit), cells = n(), votes = sum(nv), .groups="drop") %>%
  left_join(wi %>% group_by(pair) %>%
              summarise(windows = n(), anchors = sum(n_anch),
                        span_Mb = round(max(mid_Mb)-min(mid_Mb),1), .groups="drop"),
            by = "pair") %>%
  left_join(v %>% count(pair, name = "raw_votes"), by = "pair")
print(as.data.frame(pg))
cat(sprintf("\n%s  %d check(s) failed\n", ifelse(FAIL == 0, "ALL CHECKS PASSED —", "PROBLEMS —"), FAIL))
