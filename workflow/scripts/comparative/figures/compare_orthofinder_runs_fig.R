#!/usr/bin/env Rscript
# compare_orthofinder_runs_fig.R -- one-page overview of compare_orthofinder_runs.py's tables.
#   A, B: what happened to each old group of 2+ genes (identical / merged / split / reshuffled), as the share of
#         genes, overall and by the old group's size; A = like for like (old root groups vs new HOGs at the
#         shared species' common ancestor), B = old vs new Orthogroups.tsv (what the screens read)
#   C:    orthologue pairs kept between the runs, per species (Jaccard of its pairs with the other shared species)
# usage: Rscript compare_orthofinder_runs_fig.R <dir with group_fates.tsv and orthologue_pairs.tsv>
# writes <dir>/overview.png and <dir>/overview.pdf. The tables themselves are the table view.
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript compare_orthofinder_runs_fig.R <compare dir>")
dir <- args[1]
fates <- fread(file.path(dir, "group_fates.tsv"))
op_fn <- file.path(dir, "orthologue_pairs.tsv")
op <- if (file.exists(op_fn)) fread(op_fn) else NULL

# ---- palette (dataviz reference palette, light mode; validated: categorical slots 1-4, adjacent pairs)
surface <- "#fcfcfb"; ink <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"
grid_col <- "#e1e0d9"; base_col <- "#c3c2b7"
fate_levels <- c("identical", "merged", "split", "reshuffled")
fate_col <- c(identical = "#2a78d6", merged = "#eb6834", split = "#1baf7a", reshuffled = "#eda100")
label_on <- c(identical = "#ffffff", merged = "#ffffff", split = "#ffffff", reshuffled = "#0b0b0b")
bar_col <- "#2a78d6"

# ---- shares of genes by fate, overall and by old-group size
size_breaks <- c(1, 2, 5, 10, 25, 100, Inf)
size_labels <- c("2 genes", "3-5", "6-10", "11-25", "26-100", "> 100")
fates[, size := cut(n_genes, size_breaks, labels = size_labels)]
share_tab <- function(d) {
  s <- d[, .(genes = sum(n_genes)), by = fate]
  s <- merge(data.table(fate = fate_levels), s, all.x = TRUE)[is.na(genes), genes := 0]
  s[, share := 100 * genes / sum(genes)]
  setkey(s, fate)[fate_levels]
}
rows_for <- function(cmp) {
  d <- fates[comparison == cmp]
  out <- list(list(label = "all groups", n = nrow(d), s = share_tab(d)))
  for (lv in size_labels) {
    dd <- d[size == lv]
    if (nrow(dd)) out[[length(out) + 1]] <- list(label = lv, n = nrow(dd), s = share_tab(dd))
  }
  out
}

fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")

draw_fates <- function(cmp, title, subtitle) {
  rows <- rows_for(cmp)
  k <- length(rows)
  par(mar = c(3.4, 6.6, 3.9, 5.0), bg = surface, cex = 1)
  plot.new(); plot.window(xlim = c(0, 100), ylim = c(k + 0.5, 0.5), xaxs = "i")
  # bar thickness: 24 px, whatever the device size
  px <- diff(grconvertY(c(0, 24), "device", "user")); half <- abs(px) / 2
  gap <- abs(diff(grconvertX(c(0, 2), "device", "user")))   # 2 px surface gap between segments
  abline(v = seq(0, 100, 25), col = grid_col, lwd = 1)
  for (i in seq_len(k)) {
    s <- rows[[i]]$s; x0 <- 0
    for (f in fate_levels) {
      w <- s[fate == f, share]
      if (w > 0) {
        rect(x0 + if (x0 > 0) gap / 2 else 0, i - half, x0 + w - if (x0 + w < 100) gap / 2 else 0, i + half,
             col = fate_col[[f]], border = NA)
        lab <- sprintf("%.0f%%", w)
        if (strwidth(lab, cex = 0.75) + 6 * gap < w) text(x0 + w / 2, i, lab, col = label_on[[f]], cex = 0.75)
      }
      x0 <- x0 + w
    }
    mtext(rows[[i]]$label, side = 2, at = i, las = 1, line = 0.4, cex = 0.85, col = if (i == 1) ink else ink2,
          font = if (i == 1) 2 else 1)
    mtext(paste0("n = ", fmt_n(rows[[i]]$n)), side = 4, at = i, las = 1, line = 0.3, cex = 0.75, col = muted)
  }
  axis(1, at = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%"), col = base_col, col.axis = muted,
       cex.axis = 0.75, tcl = -0.25, mgp = c(2, 0.35, 0), lwd = 0, lwd.ticks = 1)
  mtext("share of genes, by what happened to their old group", side = 1, line = 2.0, cex = 0.75, col = muted)
  mtext("old group size", side = 2, line = 5.6, at = 0.25 + k / 2, cex = 0.75, col = muted)
  mtext(title, side = 3, line = 2.3, adj = 0, cex = 1.0, col = ink, font = 2)
  mtext(subtitle, side = 3, line = 1.1, adj = 0, cex = 0.8, col = ink2)
}

draw_pairs <- function() {
  sp <- unique(c(op$species_1, op$species_2))
  per <- rbindlist(lapply(sp, function(s) {
    d <- op[species_1 == s | species_2 == s]
    data.table(species = s, jac = sum(d$both) / sum(d$old_pairs + d$new_pairs - d$both))
  }))
  overall <- sum(op$both) / sum(op$old_pairs + op$new_pairs - op$both)
  setorder(per, -jac)
  k <- nrow(per)
  par(mar = c(3.4, 11.5, 3.9, 5.0), bg = surface, cex = 1)
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(k + 0.5, 0.5), xaxs = "i")
  half <- abs(diff(grconvertY(c(0, 24), "device", "user"))) / 2
  abline(v = seq(0, 1, 0.25), col = grid_col, lwd = 1)
  for (i in seq_len(k)) {
    rect(0, i - half, per$jac[i], i + half, col = bar_col, border = NA)
    text(per$jac[i], i, sprintf(" %.2f", per$jac[i]), adj = 0, cex = 0.75, col = ink2, xpd = NA)
    mtext(gsub("_", " ", per$species[i]), side = 2, at = i, las = 1, line = 0.4, cex = 0.85, col = ink2, font = 3)
  }
  axis(1, at = seq(0, 1, 0.25), col = base_col, col.axis = muted, cex.axis = 0.75, tcl = -0.25,
       mgp = c(2, 0.35, 0), lwd = 0, lwd.ticks = 1)
  mtext("Jaccard of the species' orthologue pairs in the two runs (1 = the same pairs)", side = 1, line = 2.0,
        cex = 0.75, col = muted)
  mtext("C  orthologue pairs kept between the runs, per species", side = 3, line = 2.3, adj = 0, cex = 1.0,
        col = ink, font = 2)
  mtext(sprintf("pairs with the other shared species, from the gene trees; all species pairs together: %.2f", overall),
        side = 3, line = 1.1, adj = 0, cex = 0.8, col = ink2)
}

draw_legend <- function() {
  par(mar = c(0, 0, 0, 0), bg = surface, cex = 1)
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.008, 0.74, "How adding species changed the grouping of the shared species' genes", adj = 0, cex = 1.3,
       col = ink, font = 2)
  x <- 0.008
  for (f in fate_levels) {
    rect(x, 0.17, x + 0.014, 0.37, col = fate_col[[f]], border = NA)
    text(x + 0.02, 0.27, f, adj = 0, cex = 0.85, col = ink2)
    x <- x + 0.035 + strwidth(f, cex = 0.85)
  }
  text(x + 0.01, 0.27, paste0("merged: in one new group with other genes | split: in several new groups ",
                              "holding nothing else | reshuffled: anything else"), adj = 0, cex = 0.75, col = muted)
}

render <- function() {
  layout(matrix(c(1, 1, 2, 3, 4, 4), nrow = 3, byrow = TRUE), heights = c(0.62, 3, 3.1))
  draw_legend()
  draw_fates("A", "A  like for like", "old root groups vs new HOGs at the shared species' ancestor")
  draw_fates("B", "B  old vs new Orthogroups.tsv", "what the screens read; the older root lumps families")
  if (!is.null(op) && nrow(op)) draw_pairs() else { plot.new(); text(0.5, 0.5, "no orthologue_pairs.tsv", col = muted) }
}

png(file.path(dir, "overview.png"), width = 11, height = 8.6, units = "in", res = 200, bg = surface)
render(); invisible(dev.off())
pdf_dev <- if (capabilities("cairo")) cairo_pdf else pdf
pdf_dev(file.path(dir, "overview.pdf"), width = 11, height = 8.6, bg = surface)
render(); invisible(dev.off())
cat("written:", file.path(dir, "overview.png"), "and overview.pdf\n")
