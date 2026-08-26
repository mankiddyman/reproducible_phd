#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value = TRUE); if (length(.p)) .libPaths(.p)
# BLOCK level on BOTH sides. Dionaea chromosomes are mosaics, so the phasing unit is
# the Dionaea syntenic block. Signed graph on Dionaea blocks:
#   DIFFERENT  = the two Dionaea blocks of one anchor (homeologs)
#   SAME       = two Dionaea blocks that attract the same Drosera block
# Two-colour it. Blocks on one chromosome with different colours = HE breakpoint.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ragg)
})
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype")
GSD <- "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"
MINV <- 4; MINF <- 0.70   # a Drosera block must cast >=MINV votes, >=MINF one way

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

## ---- the two Dionaea genes per anchor, and their blocks ---------------------
dio <- meta %>% filter(genome == "Dionaea_muscipula") %>%
  left_join(bed %>% transmute(gene = id, genome, mid = (start+end)/2),
            by = c("gene","genome")) %>%
  left_join(anch, by = "nep_gene") %>% filter(!is.na(mid), !is.na(nep_chr)) %>%
  group_by(nep_gene) %>% filter(n() == 2, n_distinct(chr) == 2) %>%
  arrange(chr, .by_group = TRUE) %>%
  summarise(nep_chr = nep_chr[1],
            geneA = gene[1], chrA = chr[1], midA = mid[1],
            geneB = gene[2], chrB = chr[2], midB = mid[2], .groups = "drop")
dio$dblkA <- mapply(ab, "Dionaea_muscipula", dio$chrA, dio$midA, dio$nep_chr)
dio$dblkB <- mapply(ab, "Dionaea_muscipula", dio$chrB, dio$midB, dio$nep_chr)
dio <- dio %>% filter(!is.na(dblkA), !is.na(dblkB)) %>%
  mutate(DA = paste0(sub("_sg[0-9]+_s", "-s", sub("^chr","c",chrA)), "#", dblkA),
         DB = paste0(sub("_sg[0-9]+_s", "-s", sub("^chr","c",chrB)), "#", dblkB))
cat(sprintf("anchors with both Dionaea copies in blocks: %d\n", nrow(dio)))
cat(sprintf("distinct Dionaea blocks: %d\n", n_distinct(c(dio$DA, dio$DB))))

## ---- join votes ------------------------------------------------------------
V <- v %>% inner_join(select(dio, anchor = nep_gene, DA, DB), by = "anchor") %>%
  mutate(dio_blk = ifelse(vote == "A", DA, DB))
cat(sprintf("votes with both sides block-resolved: %d\n", nrow(V)))

## ---- Drosera block -> which Dionaea block does it ally with? ---------------
R <- V %>% group_by(unit) %>%
  summarise(nv = n(), .groups = "drop_last") %>% ungroup()
RB <- V %>% count(unit, dio_blk, name = "k") %>%
  group_by(unit) %>% mutate(tot = sum(k), f = k/tot) %>% ungroup() %>%
  filter(tot >= MINV, f >= MINF)
cat(sprintf("Drosera blocks with a confident ally (>=%d votes, >=%.0f%% one way): %d\n",
            MINV, 100*MINF, n_distinct(RB$unit)))

## ---- signed graph ----------------------------------------------------------
diff_e <- dio %>% distinct(a = DA, b = DB) %>% mutate(same = FALSE)
same_e <- RB %>% group_by(unit) %>% filter(n() >= 2) %>%
  summarise(bl = list(sort(unique(dio_blk))), .groups = "drop") %>%
  rowwise() %>% do({ b <- .$bl[[1]]
    if (length(b) < 2) tibble() else {
      cb <- combn(b, 2); tibble(a = cb[1,], b = cb[2,], same = TRUE) } }) %>%
  ungroup() %>% distinct()
ed <- bind_rows(diff_e, same_e) %>% filter(a != b) %>% distinct()
cat(sprintf("edges: %d DIFFERENT, %d SAME\n", sum(!ed$same), sum(ed$same)))

nodes <- sort(unique(c(ed$a, ed$b)))
adj <- split(seq_len(nrow(ed)), NA)  # placeholder
col <- setNames(rep(NA_integer_, length(nodes)), nodes)
comp <- setNames(rep(NA_integer_, length(nodes)), nodes)
nb_of <- lapply(nodes, function(x) which(ed$a == x | ed$b == x))
names(nb_of) <- nodes
ci <- 0L
for (s in nodes) {
  if (!is.na(col[s])) next
  ci <- ci + 1L; col[s] <- 1L; comp[s] <- ci; q <- s
  while (length(q)) {
    x <- q[1]; q <- q[-1]
    for (i in nb_of[[x]]) {
      y <- if (ed$a[i] == x) ed$b[i] else ed$a[i]
      want <- if (ed$same[i]) col[x] else 3L - col[x]
      if (is.na(col[y])) { col[y] <- want; comp[y] <- ci; q <- c(q, y) }
    }
  }
}
viol <- sum(mapply(function(a,b,s) (col[a] == col[b]) != s, ed$a, ed$b, ed$same))
cat(sprintf("\n=== two-colouring: %d/%d edges violated | %d connected components ===\n",
            viol, nrow(ed), max(comp, na.rm = TRUE)))
cat(sprintf("component sizes: %s\n", paste(sort(table(comp), decreasing = TRUE), collapse = " ")))

out <- tibble(dio_block = nodes, colour = ifelse(col[nodes] == 1, "A", "B"),
              component = comp[nodes]) %>%
  mutate(chr = sub("#.*", "", dio_block), blkID = sub(".*#", "", dio_block)) %>%
  left_join(nb %>% filter(sp == "Dionaea_muscipula") %>%
              transmute(blkID = as.character(blkID), sp_chr, start = sp_s, end = sp_e),
            by = "blkID") %>% arrange(chr, start)
write_csv(out, "dionaea_block_phase.csv")

cat("\n=== do blocks on one Dionaea chromosome get different colours? ===\n")
mix <- out %>% filter(!is.na(start)) %>% group_by(chr, component) %>%
  summarise(n_blk = n(), nA = sum(colour == "A"), nB = sum(colour == "B"),
            mixed = nA > 0 & nB > 0, .groups = "drop")
print(as.data.frame(mix))
cat(sprintf("\n%d/%d chromosome-components are MIXED => HE breakpoints\n",
            sum(mix$mixed), nrow(mix)))

cat("\n=== per-chromosome block map (the coordinate dictionary) ===\n")
print(as.data.frame(out %>% filter(!is.na(start)) %>%
  transmute(chr, start_Mb = round(start/1e6,2), end_Mb = round(end/1e6,2),
            blkID, component, colour)), digits = 4)

p <- out %>% filter(!is.na(start)) %>%
  ggplot(aes(xmin = start/1e6, xmax = end/1e6,
             ymin = 0, ymax = 1, fill = colour)) +
  geom_rect(colour = "white", linewidth = .2) +
  facet_wrap(~ chr, ncol = 2, scales = "free_x") +
  scale_fill_manual(values = c(A = "#1b9e77", B = "#d95f02")) +
  labs(title = "Dionaea subgenome assignment per syntenic block",
       subtitle = "colour = inferred subgenome; a colour switch within a chromosome is an HE breakpoint",
       x = "position (Mb)", y = NULL, fill = "subgenome") +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("FIG16_dionaea_block_phase.png", p, width = 12, height = 10, dpi = 170, device = agg_png)
cat("\nWROTE: FIG16_dionaea_block_phase.png dionaea_block_phase.csv\n")
