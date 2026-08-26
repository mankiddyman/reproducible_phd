#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# GENESPACE splits one collinear tract into many small blocks. Merge blocks that are
# ADJACENT on the Drosera chromosome (gap <= MAXGAP) AND agree on which Dionaea
# homeolog they vote for. Concordance is required, so a merge cannot erase a real
# exchange boundary -- discordant neighbours stay separate.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
MAXGAP <- 5e6; MINV <- 3; MAXDIFF <- 0.34

blk <- read_csv(file.path(GSD,"syntenicBlock_coordinates.csv"), show_col_types = FALSE)
v <- read_csv("tract_votes_blocks7.csv", show_col_types = FALSE) %>%
  mutate(bid = sub(".*:\\s*", "", blk))

## block extent on the Drosera side
ext <- blk %>% filter(genome1=="Nepenthes_gracilis" | genome2=="Nepenthes_gracilis") %>%
  mutate(n1 = genome1=="Nepenthes_gracilis",
         sp = ifelse(n1,genome2,genome1), nep_chr = ifelse(n1,chr1,chr2),
         sp_chr = ifelse(n1,chr2,chr1),
         s = ifelse(n1,pmin(startBp2,endBp2),pmin(startBp1,endBp1)),
         e = ifelse(n1,pmax(startBp2,endBp2),pmax(startBp1,endBp1))) %>%
  filter(sp != "Nepenthes_gracilis") %>%
  transmute(sp, nep_chr, sp_chr, bid = as.character(blkID), s, e)

B <- v %>% group_by(pair, sp, lin_chr, bid) %>%
  summarise(nv = n(), fA = mean(vote == "A"), .groups = "drop") %>%
  left_join(ext, by = c("sp", "lin_chr" = "sp_chr", "bid")) %>%
  filter(!is.na(s)) %>% arrange(pair, sp, lin_chr, s)

## chain adjacent + concordant blocks
B <- B %>% group_by(pair, sp, lin_chr) %>%
  mutate(gap = s - lag(e),
         dif = abs(fA - lag(fA)),
         newgrp = is.na(gap) | gap > MAXGAP |
                  (lag(nv) >= MINV & nv >= MINV & !is.na(dif) & dif > MAXDIFF),
         mblk = cumsum(newgrp)) %>% ungroup()

MAP <- B %>% transmute(pair, sp, lin_chr, bid,
                       mbid = sprintf("m%02d", mblk))
write_csv(MAP, "block_merge_map.csv")

cat("=== merge summary ===\n")
print(as.data.frame(B %>% group_by(pair, sp, lin_chr) %>%
  summarise(blocks_in = n(), blocks_out = n_distinct(mblk),
            votes = sum(nv), .groups = "drop") %>%
  group_by(sp) %>% summarise(chr_tracks = n(), blocks_in = sum(blocks_in),
                             blocks_out = sum(blocks_out), votes = sum(votes),
                             .groups = "drop")))

vm <- v %>% left_join(MAP, by = c("pair","sp","lin_chr","bid")) %>%
  mutate(blk = ifelse(is.na(mbid), bid, mbid)) %>% select(-mbid, -bid)
write_csv(vm, "tract_votes_blocks7m.csv")

cat("\n=== rows that would now pass MINROW>=6 and >=2 windows ===\n")
cmp <- bind_rows(
  v  %>% group_by(pair, sp, lin_chr, bid) %>% summarise(nv = n(), .groups="drop") %>%
         mutate(set = "before") %>% rename(b = bid),
  vm %>% group_by(pair, sp, lin_chr, blk) %>% summarise(nv = n(), .groups="drop") %>%
         mutate(set = "after")  %>% rename(b = blk))
print(as.data.frame(cmp %>% group_by(set) %>%
  summarise(blocks = n(), pass6 = sum(nv >= 6), votes_in_passing = sum(nv[nv >= 6]),
            .groups = "drop")))

cat("\n=== binata on Dionaea pair chr5, after merging ===\n")
print(as.data.frame(vm %>% filter(sp=="Drosera_binata", pair=="chr5") %>%
  count(lin_chr, blk, name="votes") %>% arrange(desc(votes))))
