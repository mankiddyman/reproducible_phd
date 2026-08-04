#!/usr/bin/env Rscript
# fig_rad21_loci.R -- locus map of RAD21/SCC1 (OG0001569) members per species
suppressWarnings(suppressPackageStartupMessages({ library(data.table); library(ggplot2) }))

repo    <- normalizePath(".")
fig_dir <- file.path(repo, "results/comparative/holocentricity/figures")
loci_f  <- file.path(repo, "results/comparative/holocentricity/arabidopsis_map/og0001569_loci.tsv")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

d <- fread(loci_f)

d[, type := ifelse(toupper(type) == "HOLO", "holocentric", "monocentric")]

# collapse haplotypes to a locus label
d[, base_scaf := scaffold]
d[, base_scaf := gsub("_hap[12]$", "", base_scaf)]
d[, base_scaf := gsub("_collapsed$", "", base_scaf)]
d[species == "Drosera_roseana", base_scaf := "single locus\n(hap-resolved)"]

# short scaffold label: handle utg / tg / h1tg / h2tg prefixes
d[, scaf_lab := base_scaf]
d[grepl("tg", scaf_lab), scaf_lab := sub("^h[12]", "", scaf_lab)]   # h2tg010111l -> tg010111l
d[grepl("^utg|^tg", scaf_lab), scaf_lab := sub("l$", "", sub("0+", "", scaf_lab))]  # tg010111l -> tg10111

sp_order <- c("Drosera_paradoxa","Drosera_regia","Drosera_roseana",
              "Drosera_aliciae","Drosera_binata","Drosera_capensis",
              "Drosera_filiformis","Drosera_tokaiensis")
d[, sp := factor(gsub("Drosera_","",species), levels = rev(gsub("Drosera_","",sp_order)))]

loci_n <- d[, .(n_loci = uniqueN(base_scaf), n_members = .N), by = .(sp, type)]
setorder(loci_n, -sp); cat("=== distinct loci per species ===\n"); print(loci_n)

# x-position per locus; stack members within a locus
d[, x := as.integer(factor(base_scaf)), by = sp]
d[, n_in_locus := .N, by = .(sp, base_scaf)]
d[, member_in_locus := seq_len(.N), by = .(sp, base_scaf)]
d[, xpos := x + (member_in_locus - 1) * 0.20]

# one label per locus (alternating height)
lab_dt <- d[d[, .I[1], by = .(sp, base_scaf)]$V1]
lab_dt[, vj := ifelse(seq_len(.N) %% 2 == 0, -2.0, -1.1), by = sp]

# tandem brackets: loci with >=2 members get a connector + "xN tandem" mark
tand <- d[n_in_locus >= 2, .(x0 = min(xpos), x1 = max(xpos), n = .N,
                             ymid = as.integer(sp)[1]), by = .(sp, base_scaf, type)]
tand[, sp_y := as.integer(sp)]

col_type <- c(holocentric = "#1D9E75", monocentric = "#BA7517")

p <- ggplot(d, aes(xpos, sp, colour = type)) +
  # tandem connector: a small horizontal bracket UNDER the stacked points
  geom_segment(data = tand,
               aes(x = x0, xend = x1, y = sp_y - 0.22, yend = sp_y - 0.22, colour = type),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(data = tand, aes(x = x0, xend = x0, y = sp_y - 0.22, yend = sp_y - 0.12, colour = type),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(data = tand, aes(x = x1, xend = x1, y = sp_y - 0.22, yend = sp_y - 0.12, colour = type),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_text(data = tand, aes(x = (x0 + x1)/2, y = sp_y - 0.34, label = paste0("\u00d7", n, " tandem"), colour = type),
            inherit.aes = FALSE, size = 2.4, show.legend = FALSE) +
  geom_point(size = 4.5) +
  geom_text(data = lab_dt, aes(label = scaf_lab, vjust = vj), size = 2.5, colour = "grey30") +
  scale_colour_manual(values = col_type, name = NULL) +
  scale_x_continuous(breaks = NULL,
                     name = "distinct loci  (bracketed = multiple copies on the SAME chromosome/scaffold = tandem)",
                     expand = expansion(mult = c(0.04, 0.08))) +
  labs(title = "RAD21/SCC1 cohesin: locus map across Drosera",
       subtitle = "each point = a gene copy at its scaffold. holocentrics = single locus; monocentrics = multiple loci (different chromosomes and/or tandem)",
       y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_line(colour = "grey92"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 10.5, colour = "grey40"),
        axis.title.x = element_text(size = 9.5, colour = "grey45"))

ggsave(file.path(fig_dir, "fig_rad21_loci.pdf"), p, width = 9.5, height = 5.2)
ggsave(file.path(fig_dir, "fig_rad21_loci.png"), p, width = 9.5, height = 5.2, dpi = 200)
cat("wrote fig_rad21_loci.{pdf,png}\n")
