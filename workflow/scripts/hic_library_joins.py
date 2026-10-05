#!/usr/bin/env python3
"""hic_library_joins.py -- does every Hi-C library (sample/date) agree on how chr1 and chr2 are joined?

If the libraries came from different tissues, plants or dates, one of them carrying the
pollen's arrangement (L1 of A attached to Q of D, L2 of D to P of A) would show it here.
Streams the re-mapped pairs (hic_remap work/pairs_all.tsv: name, chrom1, pos1, mapq1,
chrom2, pos2, mapq2, library) once and counts contacts per library in 5 x 5 Mb boxes
(1 Mb left out either side of a join), as contacts per Mb x Mb, divided by the same
library's ordinary stretch of chromosome. ~1 = joined; ~background = not joined.

Usage: hic_library_joins.py PAIRS.tsv OUT.txt
"""
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

PAIRS, OUT = sys.argv[1:3]
J = lambda c, x: ((c, x - 6, x - 1), (c, x + 1, x + 6))          # boxes either side of a point on one chromosome
BOXES = {
    "A join  L1|P": ((("chr1_hap1", 257, 262), ("chr1_hap1", 264, 269)),),
    "D join  L2|Q": ((("chr2_hap2", 209, 214), ("chr2_hap2", 216, 221)),),
    "L1 of A x Q of D": ((("chr1_hap1", 257, 262), ("chr2_hap2", 216, 221)),),
    "L2 of D x P of A": ((("chr2_hap2", 209, 214), ("chr1_hap1", 264, 269)),),
    "chr5 x chr6": ((("chr5_hap1", 60, 70), ("chr6_hap1", 180, 190)),),          # the stretch the pollen link
    "control": tuple(J(c, x) for c, x in (("chr1_hap1", 100), ("chr1_hap1", 180), ("chr1_hap1", 330),
                                          ("chr2_hap2", 100), ("chr2_hap2", 270), ("chr3_hap1", 150),
                                          ("chr4_hap1", 250), ("chr6_hap1", 300))),
    "background": ((("chr1_hap1", 257, 262), ("chr5_hap1", 100, 105)), (("chr2_hap2", 209, 214), ("chr6_hap1", 200, 205)),
                   (("chr3_hap1", 144, 149), ("chr6_hap1", 300, 305)), (("chr4_hap1", 244, 249), ("chr5_hap1", 200, 205))),
}
cnt = {m: defaultdict(lambda: defaultdict(float)) for m in ("unique", "all")}   # mode -> library -> box -> contacts
area = {k: sum((a[2] - a[1]) * (b[2] - b[1]) for a, b in v) for k, v in BOXES.items()}
tot = defaultdict(int)
for ch in pd.read_csv(PAIRS, sep="\t", header=None, usecols=[1, 2, 3, 4, 5, 6, 7],
                      names=["c1", "p1", "q1", "c2", "p2", "q2", "lib"], chunksize=10_000_000,
                      dtype={"c1": str, "c2": str, "lib": str}):
    for L, n in ch.lib.value_counts().items():
        tot[L] += int(n)
    c1, c2 = ch.c1.to_numpy(), ch.c2.to_numpy()
    p1, p2 = ch.p1.to_numpy() / 1e6, ch.p2.to_numpy() / 1e6
    uq = ((ch.q1 >= 30) & (ch.q2 >= 30)).to_numpy()
    lib = ch.lib.to_numpy()
    for name, boxes in BOXES.items():
        hit = np.zeros(len(ch), dtype=bool)
        for (ca, a0, a1), (cb, b0, b1) in boxes:
            hit |= (c1 == ca) & (p1 >= a0) & (p1 < a1) & (c2 == cb) & (p2 >= b0) & (p2 < b1)
            hit |= (c2 == ca) & (p2 >= a0) & (p2 < a1) & (c1 == cb) & (p1 >= b0) & (p1 < b1)
        for mode, sel in (("all", hit), ("unique", hit & uq)):
            if sel.any():
                for L, n in pd.Series(lib[sel]).value_counts().items():
                    cnt[mode][L][name] += n
out = ["HI-C JOINS PER LIBRARY  (contacts per Mb x Mb, then / the library's ordinary stretch of chromosome)",
       "  ~1 = joined like a chromosome; ~background = not joined", ""]
keys = ["A join  L1|P", "D join  L2|Q", "L1 of A x Q of D", "L2 of D x P of A", "chr5 x chr6", "background"]
for mode in ("unique", "all"):
    out.append("%s reads (MAPQ %s both mates)" % (mode.upper(), ">= 30" if mode == "unique" else ">= 1"))
    out.append("  %-28s %11s %9s | %s" % ("library", "pairs", "control", "  ".join("%-17s" % k for k in keys)))
    for L in sorted(tot):
        c = cnt[mode][L]
        ctl = c["control"] / area["control"]
        cells = []
        for k in keys:
            v = c[k] / area[k]
            cells.append("%-17s" % ("%.3f (%d)" % (v / ctl, c[k]) if ctl > 0 else "n/a (%d)" % c[k]))
        out.append("  %-28s %11d %9.2f | %s" % (L[:28], tot[L], ctl, "  ".join(cells)))
    out.append("")
out.append("  (each cell: ratio to the control, then the raw contact count in brackets)")
open(OUT, "w").write("\n".join(out) + "\n")
print("\n".join(out))
