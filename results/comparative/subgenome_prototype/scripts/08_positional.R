#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2); library(ragg)})
BASE <- Sys.getenv("SUBG_BASE", getwd())
SET  <- commandArgs(TRUE)[1]; if (is.na(SET)) SET <- "cds"
QC   <- file.path(BASE,"trees","qc")
GS   <- if (!is.na(commandArgs(TRUE)[2])) commandArgs(TRUE)[2] else file.path(dirname(BASE), "genespace", "results", "combBed.txt")
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
roll <- function(x, w) { n <- length(x); vapply(seq_len(n), function(i)
  mean(x[max(1,i-w%/%2):min(n,i+w%/%2)]), numeric(1)) }

hr("0. combBed")
cat("  path:", GS, "\n")
if (!file.exists(GS)) stop("combBed.txt not found — pass its path as arg 2")
cb <- read_tsv(GS, show_col_types = FALSE)
cat("  columns:", paste(names(cb), collapse=", "), "\n")
idc  <- intersect(c("id","gene"), names(cb))[1]
ordc <- intersect(c("ord","order","genome_ord","chrOrd"), names(cb))[1]
cat(sprintf("  using id column '%s', order column '%s'\n", idc, ordc))
if (is.na(ordc)) { cb$ord <- NA_integer_
  cb <- cb %>% group_by(genome, chr) %>% arrange(start, .by_group=TRUE) %>%
        mutate(ord = row_number()) %>% ungroup(); ordc <- "ord"
  cat("  no order column — derived from start coordinate\n") }
pos <- cb %>% transmute(gene = .data[[idc]], genome, chr_gs = chr, ord = .data[[ordc]])

hr("1. joins")
v    <- read_csv(file.path(QC, paste0("xy_votes_", SET, ".csv")), show_col_types=FALSE)
meta <- read_tsv(file.path(BASE,"wgd7","tip_meta.tsv"), show_col_types=FALSE) %>%
        mutate(tip = gsub("@","_",tip))
v <- v %>% left_join(meta %>% select(anchor=nep_gene, tip, gene), by=c("anchor","tip"))
jr <- function(d, col) round(mean(!is.na(d[[col]])), 4)
v2 <- v %>% left_join(pos, by=c("gene","genome"))
cat(sprintf("  votes -> tip_meta gene : %.4f\n", jr(v, "gene")))
cat(sprintf("  votes -> combBed ord   : %.4f\n", jr(v2, "ord")))
if (jr(v2,"ord") < 0.9) {
  cat("  retrying with .t1 stripped\n")
  v2 <- v %>% mutate(gene2 = sub("\\.t[0-9]+$","",gene)) %>%
        left_join(pos %>% mutate(gene2 = sub("\\.t[0-9]+$","",gene)) %>% select(-gene),
                  by=c("gene2","genome"))
  cat(sprintf("  after strip            : %.4f\n", jr(v2,"ord")))
}

## Nepenthes anchor position
nep <- meta %>% filter(genome=="Nepenthes_gracilis") %>% select(anchor=nep_gene, gene) %>%
  left_join(pos %>% filter(genome=="Nepenthes_gracilis") %>% select(gene, nep_ord=ord), by="gene")
if (mean(!is.na(nep$nep_ord)) < 0.9)
  nep <- meta %>% filter(genome=="Nepenthes_gracilis") %>%
    transmute(anchor=nep_gene, gene2=sub("\\.t[0-9]+$","",gene)) %>%
    left_join(pos %>% filter(genome=="Nepenthes_gracilis") %>%
              transmute(gene2=sub("\\.t[0-9]+$","",gene), nep_ord=ord), by="gene2")
cat(sprintf("  anchors -> Nepenthes ord: %.4f\n", mean(!is.na(nep$nep_ord))))
v2 <- v2 %>% left_join(nep %>% select(anchor, nep_ord), by="anchor")

hr("2. PANEL A — orientation audit along the Nepenthes reference")
cat("  A contiguous dip present in ALL Drosera at once = Dionaea HE tract\n")
cat("  (the X/Y frame is inverted there). One species only = Drosera block.\n\n")
A <- v2 %>% filter(!is.na(nep_ord)) %>%
  group_by(region, nep_ord, anchor) %>%
  summarise(fx_all = mean(votes_X), n_tips = n(),
            n_sp = n_distinct(genome),
            n_sp_low = sum(tapply(votes_X, genome, mean) < 0.5), .groups="drop") %>%
  arrange(region, nep_ord) %>% group_by(region) %>%
  mutate(fx_roll = roll(fx_all, 15), idx = row_number()) %>% ungroup()
print(as.data.frame(A %>% group_by(region) %>%
  summarise(anchors=n(), med_fx=round(median(fx_all),3),
            min_roll=round(min(fx_roll),3), max_roll=round(max(fx_roll),3),
            below_.35 = sum(fx_roll < 0.35), .groups="drop")), row.names=FALSE)
cat("\n  candidate HE tracts (rolling fx < 0.35 AND >=4 of 5 species below 0.5):\n")
cand <- A %>% filter(fx_roll < 0.35, n_sp_low >= 4, n_sp >= 4) %>%
  group_by(region) %>% summarise(n_anchors=n(),
    ord_range = paste0(min(nep_ord), "-", max(nep_ord)), .groups="drop")
if (nrow(cand)) print(as.data.frame(cand), row.names=FALSE) else cat("    none\n")
write_csv(A, file.path(QC, paste0("panelA_orientation_", SET, ".csv")))

pA <- ggplot(A, aes(idx, fx_roll)) +
  geom_hline(yintercept=c(0.5), linetype="dashed", colour="grey50") +
  geom_point(aes(y=fx_all), alpha=0.25, size=0.7) +
  geom_line(colour="#D85A30", linewidth=0.8) +
  facet_wrap(~region, scales="free_x") + ylim(0,1) +
  labs(title="Panel A — X-vote fraction along the Nepenthes reference (all Drosera pooled)",
       subtitle="Contiguous dips shared by all species = Dionaea HE, where the X/Y frame is inverted",
       x="anchor index along Nepenthes chromosome", y="fraction voting X-side") +
  theme_minimal(10)
ggsave(file.path(QC, paste0("FIG47_panelA_", SET, ".pdf")), pA, width=11, height=7)

hr("3. PANEL B — block structure along each Drosera chromosome")
B <- v2 %>% filter(!is.na(ord)) %>%
  group_by(genome, chr_gs) %>% filter(n() >= 30) %>%
  arrange(ord, .by_group=TRUE) %>%
  mutate(fx_roll = roll(as.numeric(votes_X), 11), idx = row_number()) %>% ungroup()
cat(sprintf("  chromosomes with >=30 anchors: %d\n\n", n_distinct(paste(B$genome,B$chr_gs))))
print(as.data.frame(B %>% group_by(genome, chr_gs) %>%
  summarise(n=n(), fx=round(mean(votes_X),3),
            roll_min=round(min(fx_roll),2), roll_max=round(max(fx_roll),2),
            swing=round(max(fx_roll)-min(fx_roll),2), .groups="drop") %>%
  arrange(desc(swing)) %>% head(20)), row.names=FALSE)
cat("\n  Large swing = the chromosome is a mosaic. Flat = one ancestry along its length.\n")
write_csv(B, file.path(QC, paste0("panelB_blocks_", SET, ".csv")))

pB <- ggplot(B, aes(idx, fx_roll)) +
  geom_hline(yintercept=0.5, linetype="dashed", colour="grey50") +
  geom_line(colour="#1D9E75", linewidth=0.7) +
  facet_wrap(~paste(sub("Drosera_","D_",genome), chr_gs), scales="free_x", ncol=5) +
  ylim(0,1) +
  labs(title="Panel B — X-vote fraction along each Drosera chromosome",
       subtitle="Plateaus = blocks. Short excursions = conversion. Read only after Panel A.",
       x="anchor index along chromosome", y="fraction voting X-side") +
  theme_minimal(9)
ggsave(file.path(QC, paste0("FIG48_panelB_", SET, ".pdf")), pB, width=13, height=9)

cat(sprintf("\nWROTE: %s/{panelA_orientation,panelB_blocks}_%s.csv  FIG47 FIG48\n", QC, SET))
