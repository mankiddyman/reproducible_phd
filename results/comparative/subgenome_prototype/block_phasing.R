#!/usr/bin/env Rscript
.libPaths(grep("micromamba_envs/smk", .libPaths(), value = TRUE))
# BLOCK-level phasing. A chromosome can be an HE mosaic; a syntenic block is the
# unit that descends from one subgenome. For each block, ask whether its genes sit
# nearer Dionaea X or Dionaea Y in the gene tree. Then:
#   - within-block consistency = is the block coherent? (your validation step)
#   - positional runs test  = is there an HE tract inside the block?
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ape); library(ggplot2); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
MINGENES <- 12

bed <- read.table(file.path(GSD, "combBed.txt"), header = TRUE, sep = "\t",
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
blk <- read_csv(file.path(GSD, "syntenicBlock_coordinates.csv"), show_col_types = FALSE)
meta <- read_tsv("wgd/tip_meta.tsv", show_col_types = FALSE)
anch <- read_csv("synteny_ortho_hits.csv", show_col_types = FALSE) %>% distinct(nep_gene, nep_chr)
ref <- read_csv("anchor_drosera_fractionation.csv", show_col_types = FALSE) %>%
  filter(anchor == "Nepenthes_gracilis")
xy <- bind_rows(transmute(ref, chr = winner, side = "X"),
                transmute(ref, chr = ifelse(winner == chrA, chrB, chrA), side = "Y"))

## ---- 1. Nepenthes-anchored blocks, species-side coordinates ----------------
nb <- blk %>% filter(genome1 == "Nepenthes_gracilis" | genome2 == "Nepenthes_gracilis") %>%
  mutate(n1 = genome1 == "Nepenthes_gracilis",
         sp      = ifelse(n1, genome2, genome1),
         nep_chr = ifelse(n1, chr1, chr2),
         sp_chr  = ifelse(n1, chr2, chr1),
         sp_s    = ifelse(n1, pmin(startBp2, endBp2), pmin(startBp1, endBp1)),
         sp_e    = ifelse(n1, pmax(startBp2, endBp2), pmax(startBp1, endBp1)),
         nhit    = ifelse(n1, nHits2, nHits1)) %>%
  filter(sp != "Nepenthes_gracilis", grepl("_dom$", nep_chr), !is.na(sp_s)) %>%
  mutate(block = paste(sp, nep_chr, sp_chr, blkID, sep = "|")) %>%
  select(block, sp, nep_chr, sp_chr, sp_s, sp_e, nhit)

cat(sprintf("Nepenthes-anchored blocks: %d\n", nrow(nb)))
cat("\n=== blocks per species per ancestral region ===\n")
print(as.data.frame(nb %>% count(sp, nep_chr) %>%
  pivot_wider(names_from = nep_chr, values_from = n, values_fill = 0)))

## sanity: are the blocks of one species in one region on distinct chromosomes?
cat("\n=== block independence check (same chr = potential tandem/contiguous) ===\n")
print(as.data.frame(nb %>% group_by(sp, nep_chr) %>%
  summarise(n_blk = n(), n_chr = n_distinct(sp_chr),
            same_chr_pairs = n_blk - n_chr, .groups = "drop") %>%
  group_by(sp) %>% summarise(blocks = sum(n_blk), on_distinct_chr = sum(n_chr),
                             sharing_a_chr = sum(same_chr_pairs), .groups = "drop")))

## ---- 2. gene -> block --------------------------------------------------------
pos <- bed %>% transmute(gene = id, genome, chr, mid = (start + end) / 2, ord)
tips <- meta %>% left_join(pos, by = c("gene", "genome")) %>%
  left_join(anch, by = "nep_gene")

assign_block <- function(sp_, chr_, mid_, nchr_) {
  if (is.na(mid_)) return(NA_character_)
  c <- nb[nb$sp == sp_ & nb$sp_chr == chr_ & nb$nep_chr == nchr_ &
          nb$sp_s <= mid_ & nb$sp_e >= mid_, ]
  if (!nrow(c)) return(NA_character_)
  c$block[which.max(c$nhit)]
}
tips$block <- NA_character_
i <- which(tips$genome != "Nepenthes_gracilis")
tips$block[i] <- mapply(assign_block, tips$genome[i], tips$chr[i],
                        tips$mid[i], tips$nep_chr[i])
cat(sprintf("\ntips assigned to a block: %.1f%% (%d of %d non-Nepenthes)\n",
            100 * mean(!is.na(tips$block[i])), sum(!is.na(tips$block[i])), length(i)))

tips <- tips %>% left_join(xy, by = "chr")

## ---- 3. per gene tree: is each tip nearer Dio X or Dio Y? -------------------
fs <- list.files("wgd/tre", pattern = "\\.tre$", full.names = TRUE)
one <- function(f) {
  tr <- tryCatch(read.tree(f), error = function(e) NULL)
  if (is.null(tr) || Ntip(tr) < 4) return(NULL)
  a <- sub("\\.tre$", "", basename(f))
  m <- tips[match(tr$tip.label, tips$tip), ]
  if (any(is.na(m$genome))) return(NULL)
  nep <- tr$tip.label[m$genome == "Nepenthes_gracilis"]
  if (length(nep) != 1) return(NULL)
  tr <- tryCatch(root(tr, nep, resolve.root = TRUE), error = function(e) NULL)
  if (is.null(tr)) return(NULL)
  m <- tips[match(tr$tip.label, tips$tip), ]
  dx <- which(m$genome == "Dionaea_muscipula" & m$side == "X")
  dy <- which(m$genome == "Dionaea_muscipula" & m$side == "Y")
  if (length(dx) != 1 || length(dy) != 1) return(NULL)
  nd <- node.depth(tr)
  j <- which(!m$genome %in% c("Nepenthes_gracilis", "Dionaea_muscipula") & !is.na(m$block))
  if (!length(j)) return(NULL)
  do.call(rbind, lapply(j, function(k) {
    ax <- nd[getMRCA(tr, tr$tip.label[c(k, dx)])]
    ay <- nd[getMRCA(tr, tr$tip.label[c(k, dy)])]
    if (ax == ay) return(NULL)
    data.frame(anchor = a, nep_chr = m$nep_chr[k], sp = m$genome[k],
               block = m$block[k], gene = m$gene[k], ord = m$ord[k],
               near = ifelse(ax < ay, "X", "Y"), stringsAsFactors = FALSE)
  }))
}
v <- bind_rows(lapply(fs, one))
write_csv(v, "block_votes.csv")
cat(sprintf("block-assigned votes: %d over %d blocks\n", nrow(v), n_distinct(v$block)))

## ---- 4. per block: assignment, consistency, HE tract test ------------------
runs_p <- function(s) {
  s <- s[!is.na(s)]; n1 <- sum(s == "X"); n2 <- sum(s == "Y"); n <- n1 + n2
  if (n1 < 4 || n2 < 4) return(NA_real_)
  r <- 1 + sum(s[-1] != s[-n])
  mu <- 2*n1*n2/n + 1; vv <- 2*n1*n2*(2*n1*n2 - n)/(n^2*(n-1))
  2 * pnorm(-abs((r - mu)/sqrt(vv)))
}

bs <- v %>% arrange(block, ord) %>% group_by(sp, nep_chr, block) %>%
  summarise(n = n(), X = sum(near == "X"), frac_X = mean(near == "X"),
            runs_p = runs_p(near), .groups = "drop") %>%
  filter(n >= MINGENES) %>%
  rowwise() %>% mutate(p = binom.test(X, n, 0.5)$p.value) %>% ungroup() %>%
  mutate(p_adj = p.adjust(p, "BH"), runs_adj = p.adjust(runs_p, "BH"),
         call = case_when(p_adj < 0.05 & frac_X > 0.5 ~ "X",
                          p_adj < 0.05 & frac_X < 0.5 ~ "Y", TRUE ~ "ambiguous"),
         chr = sub("^[^|]*\\|[^|]*\\|([^|]*)\\|.*$", "\\1", block)) %>%
  arrange(sp, nep_chr, desc(frac_X))
write_csv(bs, "block_subgenome_assignment.csv")

cat("\n=== block assignments ===\n")
print(as.data.frame(bs %>% select(sp, nep_chr, chr, n, frac_X, p_adj, call, runs_adj)),
      digits = 3)
cat("\n=== calls per species ===\n")
print(as.data.frame(count(bs, sp, call) %>%
  pivot_wider(names_from = call, values_from = n, values_fill = 0)))
cat(sprintf("\nblocks with a significant HE tract (runs_adj<0.05): %d of %d\n",
            sum(bs$runs_adj < 0.05, na.rm = TRUE), sum(!is.na(bs$runs_adj))))

## do blocks on ONE chromosome agree? (chromosome coherence, measured not assumed)
cat("\n=== chromosomes with >1 block: do the blocks agree? ===\n")
mc <- bs %>% group_by(sp, chr) %>% filter(n() > 1) %>%
  summarise(n_blk = n(), calls = paste(sort(unique(call)), collapse = "/"),
            spread = round(max(frac_X) - min(frac_X), 3), .groups = "drop")
print(as.data.frame(mc), digits = 3)

## ---- 5. plots ---------------------------------------------------------------
p1 <- ggplot(bs, aes(frac_X, reorder(paste0(chr, " (", nep_chr, ")"), frac_X),
                     colour = call, size = n)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey40") +
  geom_point() + facet_wrap(~ sp, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c(X = "#1b9e77", Y = "#d95f02", ambiguous = "grey60")) +
  scale_size_continuous(range = c(1.5, 4)) +
  labs(title = "Block-level subgenome assignment",
       subtitle = "each point = one syntenic block; X/Y = which Dionaea subgenome it sits nearer",
       x = "fraction of genes nearer Dionaea X", y = NULL, colour = "call", size = "genes") +
  theme_bw(base_size = 10)
ggsave("FIG4_block_assignment.png", p1, width = 10, height = 11, dpi = 150, device = agg_png)
ggsave("FIG4_block_assignment.pdf", p1, width = 10, height = 11)

top <- bs %>% group_by(sp) %>% slice_max(n, n = 8) %>% ungroup()
p2 <- v %>% semi_join(top, by = "block") %>% arrange(block, ord) %>%
  group_by(block) %>% mutate(idx = row_number(),
    lab = paste0(sub("^([^|]*)\\|.*$", "\\1", block[1]), " ",
                 sub("^[^|]*\\|[^|]*\\|([^|]*)\\|.*$", "\\1", block[1]))) %>% ungroup() %>%
  ggplot(aes(idx, 1, fill = near)) + geom_tile() +
  facet_wrap(~ lab, scales = "free_x", ncol = 4) +
  scale_fill_manual(values = c(X = "#1b9e77", Y = "#d95f02")) +
  labs(title = "Gene-by-gene assignment along each block",
       subtitle = "contiguous runs of one colour = homoeologous exchange tract; salt-and-pepper = noise",
       x = "gene order within block", y = NULL, fill = NULL) +
  theme_bw(base_size = 9) + theme(axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank())
ggsave("FIG5_block_positional.png", p2, width = 13, height = 7, dpi = 150, device = agg_png)
ggsave("FIG5_block_positional.pdf", p2, width = 13, height = 7)
cat("\nWROTE: FIG4_block_assignment.{png,pdf} FIG5_block_positional.{png,pdf}\n")
cat("WROTE: block_votes.csv block_subgenome_assignment.csv\n")
