#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# Map every gene-tree tip to its syntenic block. Blocks become supermatrix taxa.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
MING <- 8

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
blk <- read_csv(file.path(GSD, "syntenicBlock_coordinates.csv"), show_col_types = FALSE)
meta <- read_tsv("wgd/tip_meta.tsv", show_col_types = FALSE)
anch <- read_csv("synteny_ortho_hits.csv", show_col_types = FALSE) %>% distinct(nep_gene, nep_chr)

nb <- blk %>% filter(genome1 == "Nepenthes_gracilis" | genome2 == "Nepenthes_gracilis") %>%
  mutate(n1 = genome1 == "Nepenthes_gracilis",
         sp = ifelse(n1, genome2, genome1), nep_chr = ifelse(n1, chr1, chr2),
         sp_chr = ifelse(n1, chr2, chr1),
         sp_s = ifelse(n1, pmin(startBp2, endBp2), pmin(startBp1, endBp1)),
         sp_e = ifelse(n1, pmax(startBp2, endBp2), pmax(startBp1, endBp1)),
         nhit = ifelse(n1, nHits2, nHits1)) %>%
  filter(sp != "Nepenthes_gracilis", grepl("_dom$", nep_chr), !is.na(sp_s)) %>%
  mutate(abbr = case_when(grepl("Dionaea", sp) ~ "Dio",
                          grepl("regia", sp) ~ "reg",
                          grepl("binata", sp) ~ "bin", TRUE ~ "oth"),
         block = paste0(abbr, "_", sub("_hap1$|_collapsed$", "",
                        sub("^chr", "c", sp_chr)), "_b", blkID)) %>%
  select(block, sp, nep_chr, sp_chr, sp_s, sp_e, nhit)

pos <- bed %>% transmute(gene = id, genome, gmid = (start + end) / 2)
tips <- meta %>% left_join(pos, by = c("gene", "genome")) %>% left_join(anch, by = "nep_gene")

ab <- function(sp_, chr_, mid_, nchr_) {
  if (is.na(mid_) || is.na(nchr_)) return(NA_character_)
  c <- nb[nb$sp == sp_ & nb$sp_chr == chr_ & nb$nep_chr == nchr_ &
          nb$sp_s <= mid_ & nb$sp_e >= mid_, ]
  if (!nrow(c)) return(NA_character_)
  c$block[which.max(c$nhit)]
}
i <- which(tips$genome != "Nepenthes_gracilis")
tips$block <- NA_character_
tips$block[i] <- mapply(ab, tips$genome[i], tips$chr[i], tips$gmid[i], tips$nep_chr[i])
tips$block[tips$genome == "Nepenthes_gracilis"] <- "NEP"

keep <- tips %>% filter(!is.na(block), block != "NEP") %>%
  count(block, name = "ngenes") %>% filter(ngenes >= MING)
cat(sprintf("blocks with >=%d genes: %d\n", MING, nrow(keep)))
print(as.data.frame(keep %>% mutate(sp = sub("_.*", "", block)) %>%
  count(sp, name = "blocks")))

out <- tips %>% filter(block %in% c("NEP", keep$block)) %>%
  select(nep_gene, nep_chr, tip, gene, genome, chr, block)
write_tsv(out, "super/tip_block.tsv")
cat(sprintf("\ntips retained: %d over %d anchors, %d regions\n",
            nrow(out), n_distinct(out$nep_gene), n_distinct(out$nep_chr)))
print(as.data.frame(out %>% filter(block != "NEP") %>%
  distinct(nep_chr, block) %>% count(nep_chr, name = "blocks")))
