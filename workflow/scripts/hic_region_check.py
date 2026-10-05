#!/usr/bin/env python3
"""hic_region_check.py -- are two pieces physically joined in the Hi-C tissue? Region-level numbers.

For each test: mean contacts per bin pair between region X and region Y, divided by
  (1) the same comparison inside one ordinary chromosome at a similar spacing  -> ~1 = joined like a chromosome
  (2) unrelated chromosomes                                                 -> ~1 = no more than background
Reads the matrices written by hic_remap_verify.py (contacts_unique.npz: haplotype-specific).

Usage: hic_region_check.py CONTACTS.npz REF.fa.fai
"""
import sys

import numpy as np

d = np.load(sys.argv[1])
M, B = d["M"], int(d["bin"])
lens = {}
for l in open(sys.argv[2]):
    f = l.split("\t")
    lens[f[0]] = int(f[1])
chroms = [str(c) for c in d["chroms"]]
off, o = {}, 0
for c in chroms:
    off[c] = o
    o += -(-lens[c] // B)


def mean(c1, a1, b1, c2, a2, b2):
    i = np.arange(off[c1] + int(a1 * 1e6 // B), off[c1] + int(b1 * 1e6 // B))
    j = np.arange(off[c2] + int(a2 * 1e6 // B), off[c2] + int(b2 * 1e6 // B))
    return float(M[np.ix_(i, j)].mean())


BG = np.median([mean("chr1_hap1", 100, 140, "chr5_hap1", 100, 140), mean("chr3_hap1", 100, 140, "chr6_hap1", 100, 140),
                mean("chr2_hap2", 50, 90, "chr4_hap1", 200, 240), mean("chr5_hap1", 20, 60, "chr3_hap1", 200, 240)])
TESTS = [  # label, region X, region Y, same-spacing control inside one chromosome
    ("chr4 piece 0-85 (hap1 copy) vs the rest of chr4_hap1", ("chr4_hap1", 30, 70), ("chr4_hap1", 95, 135),
     (("chr4_hap1", 150, 190), ("chr4_hap1", 215, 255))),
    ("chr4 piece 0-85 (hap1 copy) vs the end of chr3_hap1", ("chr4_hap1", 30, 70), ("chr3_hap1", 230, 270),
     (("chr4_hap1", 150, 190), ("chr4_hap1", 215, 255))),
    ("chr5_hap1 middle vs chr6_hap1 middle", ("chr5_hap1", 100, 160), ("chr6_hap1", 150, 210),
     (("chr5_hap1", 20, 80), ("chr5_hap1", 160, 220))),
    ("chr5_hap2 middle vs chr6_hap2 middle", ("chr5_hap2", 100, 160), ("chr6_hap2", 150, 210),
     (("chr5_hap2", 20, 80), ("chr5_hap2", 160, 220))),
    ("chr5_hap1 middle vs chr6_hap2 middle", ("chr5_hap1", 100, 160), ("chr6_hap2", 150, 210),
     (("chr5_hap1", 20, 80), ("chr5_hap1", 160, 220))),
    ("chr5_hap2 middle vs chr6_hap1 middle", ("chr5_hap2", 100, 160), ("chr6_hap1", 150, 210),
     (("chr5_hap2", 20, 80), ("chr5_hap2", 160, 220))),
]
print("HI-C REGION CHECK  (%s; %d kb bins)  unrelated-chromosome background %.4f per bin pair" % (sys.argv[1], B / 1e3, BG))
print("  %-52s %10s %18s %16s" % ("test", "per bin pair", "x same-spacing cis", "x background"))
for lab, x, y, (cx, cy) in TESTS:
    if not all(r[0] in off for r in (x, y, cx, cy)):
        continue
    v, c = mean(*x, *y), mean(*cx, *cy)
    print("  %-52s %10.4f %18.2f %16.1f" % (lab, v, v / c if c else float("nan"), v / BG if BG else float("nan")))
print("  (~1 x same-spacing cis = joined like one chromosome; ~1 x background = not joined)")
