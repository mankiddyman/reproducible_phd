#!/usr/bin/env python3
"""haphic_links_check.py -- which final chromosomes does the full Hi-C data tie together?

Uses the Hi-C link counts HapHiC computed from ALL the Hi-C reads it was given
(pass-2 01.cluster/full_links.pkl: links between its input units, the RagTag units),
maps every unit to the final chromosome it ended up in (final AGP, chromosome
names by length from the .fai), and sums the links between every pair of final
chromosomes, per Mb x Mb.

Reading it: two scaffolds whose pieces really sit on one chromosome (a misjoin)
light up like a chromosome; scaffolds sharing allelic copies of a piece (hap1 vs
hap2 of chr1, or chr1_hap1 vs chr2_hap1, which both carry P) are warm from reads
that fit either copy; everything else is background.

Usage: haphic_links_check.py FULL_LINKS.pkl RAGTAG.agp FINAL.agp REF.fa.fai OUTDIR [HAPHIC_SCRIPTS_DIR]
Writes OUTDIR/haphic_links_chrom.png and OUTDIR/haphic_links_chrom.txt
"""
import os
import pickle
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

pkl, rag_agp, fin_agp, fai, out = sys.argv[1:6]
if len(sys.argv) > 6:
    sys.path.insert(0, sys.argv[6])          # in case the pickle refers to HapHiC's own classes
os.makedirs(out, exist_ok=True)

# which pieces of the chr1/chr2 scaffolds are allelic copies of each other (shared arms)
SHARED = {frozenset(("chr1_hap1", "chr1_hap2")): "both carry L1", frozenset(("chr2_hap1", "chr2_hap2")): "both carry L2",
          frozenset(("chr1_hap1", "chr2_hap1")): "both carry P", frozenset(("chr1_hap2", "chr2_hap2")): "both carry Q",
          frozenset(("chr3_hap2", "chr4_hap1")): "both carry chr4 0-85 Mb"}
DECISIVE = frozenset(("chr1_hap1", "chr2_hap2"))   # hot only if L1 sits with Q and L2 with P

links = pickle.load(open(pkl, "rb"))
items = list(links.items()) if hasattr(links, "items") else []
print("full_links.pkl: %s, %d entries; first three: %s" % (type(links).__name__, len(items), items[:3]))
if not items or not isinstance(items[0][0], tuple) or not all(isinstance(x, str) for x in items[0][0]):
    sys.exit("unexpected format: keys are not (unit, unit) name pairs -- paste the line above")

unit_len = defaultdict(int)
for l in open(rag_agp):
    if l.startswith("#") or not l.strip():
        continue
    f = l.split("\t")
    unit_len[f[0]] = max(unit_len[f[0]], int(f[2]))
chrom_by_len = {}
for l in open(fai):
    f = l.split("\t")
    chrom_by_len[int(f[1])] = f[0]
fin = defaultdict(list)
for l in open(fin_agp):
    if l.startswith("#") or not l.strip():
        continue
    f = l.rstrip("\n").split("\t")
    fin[f[0]].append(f)
obj_chrom = {s: chrom_by_len.get(max(int(f[2]) for f in fs), "unplaced") for s, fs in fin.items()}
where = defaultdict(lambda: defaultdict(int))           # unit -> chromosome -> bp placed there
span = {}                                               # unit -> (chromosome, start, end) of its largest piece
for s, fs in fin.items():
    for f in fs:
        if f[4] in ("N", "U"):
            continue
        c, n = obj_chrom[s], int(f[2]) - int(f[1]) + 1
        where[f[5]][c] += n
        if f[5] not in span or n > span[f[5]][3]:
            span[f[5]] = (c, int(f[1]), int(f[2]), n)
unit_chrom = {u: max(d, key=d.get) for u, d in where.items()}
split = [u for u, d in where.items() if len(d) > 1]
if split:
    print("units split over several chromosomes (assigned to the largest piece): %s" % ", ".join(split[:10]))

chroms = sorted({c for c in unit_chrom.values() if c != "unplaced"},
                key=lambda c: (int(c.split("_")[0][3:]), c))
idx = {c: i for i, c in enumerate(chroms)}
L = np.zeros(len(chroms))
for u, c in unit_chrom.items():
    if c in idx:
        L[idx[c]] += unit_len.get(u, where[u][c]) / 1e6
M = np.zeros((len(chroms), len(chroms)))
unknown = 0
for (u, v), n in items:
    cu, cv = unit_chrom.get(u), unit_chrom.get(v)
    if cu is None or cv is None:
        unknown += 1
        continue
    if cu in idx and cv in idx and u != v:
        M[idx[cu], idx[cv]] += n
        if idx[cu] != idx[cv]:
            M[idx[cv], idx[cu]] += n
D = M / np.outer(L, L)                                  # links per Mb x Mb
np.fill_diagonal(D, np.nan)
print("unit pairs whose units are not in the final AGP: %d" % unknown)

off = [(D[i, j], chroms[i], chroms[j]) for i in range(len(chroms)) for j in range(i + 1, len(chroms))]
vals = np.array([v for v, _, _ in off])
bg = np.nanmedian(vals)
rows = ["HI-C LINKS BETWEEN FINAL CHROMOSOMES  (HapHiC full_links, all Hi-C reads; links per Mb x Mb)",
        "  background = median over all %d chromosome pairs: %.3g" % (len(off), bg), ""]
for v, a, b in sorted(off, reverse=True)[:16]:
    tag = SHARED.get(frozenset((a, b)), "")
    if frozenset((a, b)) == DECISIVE:
        tag = "<-- DECISIVE: hot only if L1 sits with Q and L2 with P (misjoin)"
    elif a.split("_")[0] == b.split("_")[0] and not tag:
        tag = "homologs"
    elif {a.split("_")[0], b.split("_")[0]} == {"chr5", "chr6"}:
        tag = "chr5/chr6: the pollen say a translocation here"
    rows.append("  %-10s %-10s %10.3g   %6.1f x background   %s" % (a, b, v, v / bg, tag))
i, j = idx.get("chr1_hap1"), idx.get("chr2_hap2")
if i is not None and j is not None:
    rank = 1 + sum(1 for v, _, _ in off if v > D[i, j])
    rows += ["", "  chr1_hap1 x chr2_hap2 (the decisive pair): %.3g = %.1f x background, rank %d of %d" % (
        D[i, j], D[i, j] / bg, rank, len(off))]
    k, m = idx.get("chr1_hap2"), idx.get("chr2_hap1")
    if k is not None and m is not None:
        rows.append("  chr1_hap2 x chr2_hap1 (its control: no shared pieces, no misjoin claimed): %.3g = %.1f x background" % (
            D[k, m], D[k, m] / bg))
open(os.path.join(out, "haphic_links_chrom.txt"), "w").write("\n".join(rows) + "\n")
print("\n".join(rows))

fig, ax = plt.subplots(figsize=(8.6, 7.6))
R = np.log10(D / bg)
im = ax.imshow(R, cmap="Reds", vmin=0, vmax=max(1.0, np.nanpercentile(R, 99)), interpolation="nearest")
ax.set_xticks(range(len(chroms)))
ax.set_xticklabels(chroms, rotation=90, fontsize=9)
ax.set_yticks(range(len(chroms)))
ax.set_yticklabels(chroms, fontsize=9)
if i is not None and j is not None:
    for a, b in ((i, j), (j, i)):
        ax.add_patch(plt.Rectangle((b - 0.5, a - 0.5), 1, 1, fill=False, ec="#0891B2", lw=2.2))
cb = fig.colorbar(im, ax=ax, shrink=0.8)
cb.set_label("log10( links per Mb x Mb / background )")
ax.set_title("Hi-C links between final chromosomes (HapHiC, all reads)", fontsize=10)
ax.text(0, -0.2, "Boxed: chr1_hap1 x chr2_hap2. Hot = their arms sit on one chromosome (misjoin); background = a real "
        "translocation.\nWarm, expected: hap1 vs hap2 of a chromosome, and scaffolds carrying copies of the same arm "
        "(chr1_hap1/chr2_hap1 share P, chr1_hap2/chr2_hap2 share Q).",
        transform=ax.transAxes, fontsize=8, color="#334155", va="top")
fig.tight_layout()
fig.savefig(os.path.join(out, "haphic_links_chrom.png"), dpi=150, bbox_inches="tight")
print("wrote %s/haphic_links_chrom.png and .txt" % os.path.abspath(out))
