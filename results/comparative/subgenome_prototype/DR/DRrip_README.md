# Riparian subgenome figure — 2026-09-17

`DRrip6_riparian_FINAL.R` is the only script to run. The rest are the
diagnostics that established its design; keep them for the methods section.

| script | role |
|---|---|
| DRrip1_probe_affine.R    | proved plot x is affine in gene order, per chromosome |
| DRrip2_subset_search.R   | failed subset search (GENESPACE plots a subset of combBed) — kept as the negative result |
| DRrip3_offset_map.R      | the answer: `plot_ord = rank(ord) - offset(chr)`, constant offset 0..98, no flips |
| DRrip4_rank_vs_bp.R      | capensis chr15 in both coordinate systems; gene-density skew |
| DRrip5_pixel_forensic.R  | palette collision + what is actually rendered at the bar edge |
| DRrip6_riparian_FINAL.R  | **production figure** |

## Coordinate map
    plot_ord = rank(ord) within chromosome - offset(chromosome)
    x        = chromosome$x1 + plot_ord
Offset is constant per chromosome (0..98), from 4184 anchor genes, zero spread.
No chromosome is flipped. Verified: 0.000 plot units of error across 2556
block sides, by gene identity against GENESPACE's own polygon attachments.

## Three defects fixed (in order of discovery)
1. **Coordinate** — `mb2x` used rank without the per-chromosome offset, and two
   other places rescaled rank by combBed gene count over box width. Three
   different maps, all wrong. Now one map, `BED$px`, used everywhere.
2. **Audit** — the old check filtered bars by ancestral region and tested only
   the braid-end midpoint, so a wrong-coloured bar under a different region was
   scored as no mismatch. Now region-agnostic and span-weighted.
3. **Visual** — bars were drawn at y-0.05..y-0.11 while braids attach at the box
   edge (capensis: 6.9385, inside that band), so the attachment was painted over
   and the first visible pixels were 0.048 lower, after sideways drift. Tracts
   are now strips at the box edges, and each ribbon leaves vertically (STEM).

## Knobs in DRrip6
    FR   0.30   strip height, fraction of chromosome box
    HAIR 0.012  white hairline between strip and braid
    STEM 0.12   vertical riser at each attachment (ribbon spans are exactly 1.0)

## Superseded
DR16_riparian.R and DR18_riparian_full.R predate this work and carry defect 1.

## Open, not chased
- 19 braids dropped for attaching across an A/B boundary (`out/DR25_dropped.csv`)
  — may be genuine chimeric blocks rather than plotting casualties.
- 424 blocks classed discordant (ends disagree on subgenome).
- 4 of 8 region hues are nearer the green A strip than the orange B strip in RGB
  distance (`out/DR23_palette.csv`), chr6_dom and chr7_dom among them. Hue and
  subgenome coding are independent; the subtitle says so.
