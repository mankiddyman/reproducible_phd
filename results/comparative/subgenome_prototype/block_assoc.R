#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# For each Dionaea BLOCK: which Drosera lineage-blocks ally with it, and how consistently?
# attraction(drosera block, dionaea block) = votes for that Dionaea block
#   / votes cast by that Drosera block across anchors where that Dionaea block was
#     one of the two homeologs.
# 1.0 = always allies with it; 0.0 = always allies with its homeolog; 0.5 = no preference.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
MINV <- 4

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
blk <- read_csv(file.path(GSD, "syntenicBlock_coordinates.csv"), show_col_types = FALSE)
anch <- read_csv("synteny_ortho_hits.csv", show_col_types = FALSE) %>% distinct(nep_gene, nep_chr)
meta <- read_tsv("wgd/tip_meta.tsv", show_col_types = FALSE)
v <- read_csv("tract_votes_blocks.csv", show_col_types = FALSE)

nb <- blk %>% filter(genome1 == "Nepenthes_gracilis" | genome2 == "Nepenthes_gracilis") %>%
  mutate(n1 = genome1 == "Nepenthes_gracilis",
         sp = ifelse(n1, genome2, genome1), nep_chr = ifelse(n1, chr1, chr2),
         sp_chr = ifelse(n1, chr2, chr1),
         sp_s = ifelse(n1, pmin(startBp2, endBp2), pmin(startBp1, endBp1)),
         sp_e = ifelse(n1, pmax(startBp2, endBp2), pmax(startBp1, endBp1)),
         nhit = ifelse(n1, nHits2, nHits1)) %>%
  filter(sp != "Nepenthes_gracilis", grepl("_dom$", nep_chr), !is.na(sp_s)) %>%
  select(blkID, sp, nep_chr, sp_chr, sp_s, sp_e, nhit)

ab <- function(sp_, chr_, mid_, nchr_) {
  if (is.na(mid_) || is.na(nchr_)) return(NA_character_)
  c <- nb[nb$sp == sp_ & nb$sp_chr == chr_ & nb$nep_chr == nchr_ &
          nb$sp_s <= mid_ & nb$sp_e >= mid_, ]
  if (!nrow(c)) return(NA_character_)
  as.character(c$blkID[which.max(c$nhit)])
}

dio <- meta %>% filter(genome == "Dionaea_muscipula") %>%
  left_join(bed %>% transmute(gene = id, genome, mid = (start+end)/2),
            by = c("gene","genome")) %>%
  left_join(anch, by = "nep_gene") %>% filter(!is.na(mid), !is.na(nep_chr)) %>%
  group_by(nep_gene) %>% filter(n() == 2, n_distinct(chr) == 2) %>%
  arrange(chr, .by_group = TRUE) %>%
  summarise(nep_chr = nep_chr[1], chrA = chr[1], midA = mid[1],
            chrB = chr[2], midB = mid[2], .groups = "drop")
dio$bA <- mapply(ab, "Dionaea_muscipula", dio$chrA, dio$midA, dio$nep_chr)
dio$bB <- mapply(ab, "Dionaea_muscipula", dio$chrB, dio$midB, dio$nep_chr)
dio <- dio %>% filter(!is.na(bA), !is.na(bB)) %>%
  mutate(DA = paste0(sub("_sg[0-9]+_s", "-s", sub("^chr","c",chrA)), "#", bA),
         DB = paste0(sub("_sg[0-9]+_s", "-s", sub("^chr","c",chrB)), "#", bB),
         dpair = paste(DA, DB, sep = " | "))

V <- v %>% select(-any_of("nep_chr")) %>%
  inner_join(select(dio, anchor = nep_gene, nep_chr, DA, DB, dpair),
             by = "anchor") %>%
  mutate(chosen = ifelse(vote == "A", DA, DB),
         dros = paste0(sub("Drosera_","",sp), " ",
                       sub("_hap1$|_collapsed$","", lin_chr), "#", blk))
cat(sprintf("votes usable: %d | Dionaea block pairs: %d | Drosera blocks: %d\n",
            nrow(V), n_distinct(V$dpair), n_distinct(V$dros)))

## ---- attraction per (Drosera block, Dionaea block pair) --------------------
A <- V %>% group_by(nep_chr, dpair, DA, DB, dros, sp) %>%
  summarise(n = n(), toA = sum(chosen == DA), att = toA/n, .groups = "drop") %>%
  filter(n >= MINV) %>% rowwise() %>%
  mutate(p = binom.test(toA, n, 0.5)$p.value) %>% ungroup() %>%
  mutate(p_adj = p.adjust(p, "BH"),
         call = case_when(att >= 0.75 ~ "allies DA", att <= 0.25 ~ "allies DB",
                          TRUE ~ "mixed"))
write_csv(A, "block_association.csv")

cat(sprintf("\n(Drosera block, Dionaea pair) associations with >=%d votes: %d\n", MINV, nrow(A)))
cat("\n=== how decisive are they? ===\n")
print(as.data.frame(count(A, call)))
cat(sprintf("median votes per association: %.0f\n", median(A$n)))

cat("\n=== per Dionaea block pair: do Drosera blocks split both ways? ===\n")
S <- A %>% group_by(nep_chr, dpair) %>%
  summarise(n_dros = n(), nDA = sum(call == "allies DA"),
            nDB = sum(call == "allies DB"), nmix = sum(call == "mixed"),
            both_sides = nDA > 0 & nDB > 0, .groups = "drop") %>%
  arrange(desc(nDA + nDB))
print(as.data.frame(head(S, 30)))
cat(sprintf("\n%d/%d Dionaea block pairs have Drosera blocks allying BOTH ways => locally phased\n",
            sum(S$both_sides), nrow(S)))

## ---- do two Drosera blocks on the SAME chromosome ally oppositely? ---------
cat("\n=== same Drosera chromosome, different blocks, same Dionaea pair ===\n")
M <- A %>% mutate(dchr = sub("#.*", "", dros)) %>%
  group_by(dpair, dchr) %>% filter(n_distinct(dros) >= 2) %>%
  summarise(n_blk = n_distinct(dros), min_att = min(att), max_att = max(att),
            spread = max_att - min_att, opposite = min_att <= .25 & max_att >= .75,
            .groups = "drop") %>% arrange(desc(spread))
print(as.data.frame(head(M, 25)), digits = 3)
if (nrow(M)) cat(sprintf("\n%d/%d have blocks allying OPPOSITE ways on one chromosome\n",
                         sum(M$opposite), nrow(M)))

## ---- heatmap, faceted by Nepenthes region ---------------------------------
p <- ggplot(A, aes(dpair, dros, fill = att)) +
  geom_tile(colour = "white", linewidth = .3) +
  geom_text(aes(label = n), size = 1.9, colour = "grey20") +
  facet_wrap(~ nep_chr, scales = "free", ncol = 2) +
  scale_fill_gradient2(low = "#d95f02", mid = "grey92", high = "#1b9e77",
                       midpoint = .5, limits = c(0, 1)) +
  labs(title = "Which Dionaea block does each Drosera block ally with?",
       subtitle = "green = allies the first Dionaea block of the pair, orange = the second; number = votes",
       x = "Dionaea homeolog block pair", y = "Drosera lineage block", fill = "attraction") +
  theme_bw(base_size = 7) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 5),
        axis.text.y = element_text(size = 5))
ggsave("FIG17_block_association.png", p, width = 15, height = 18, dpi = 160, device = agg_png)
cat("\nWROTE: FIG17_block_association.png block_association.csv\n")
