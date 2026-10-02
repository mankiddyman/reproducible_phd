#!/usr/bin/env python3
"""hic_remap_verify.py -- every paradoxa Hi-C read re-mapped onto the final 12 chromosomes:
is each library usable, and do the contacts support each suspect join?

Input, one line per read pair (mates mapped on their own with minimap2 -x sr, primary hit,
MAPQ >= 1): name, chrom1, pos1, mapq1, chrom2, pos2, mapq2, library.
Mates < 2 kb apart on one chromosome are not Hi-C ligations and are dropped. Two contact
matrices: 'all' (MAPQ >= 1 both mates) and 'unique' (MAPQ >= 30 both: haplotype-specific).

Join score = contacts between the 5 Mb either side of a join (1 Mb left out each side) /
the same box inside ordinary chromosome arms. ~1: joined in the nucleus; near the
unrelated-chromosome level: not joined. Boxes expected to hold < 50 contacts are not scored.

Writes to OUT: library_stats.txt, join_scores.txt, hicmap_chr1_chr2.png, hicmap_chr1_chr2_zoom.png,
hicmap_chr1_chr2_zoom_unique.png, hicmap_chr1_chr2_zoom_control.png, hicmap_chr3_chr4.png,
hicmap_chr5_chr6.png, bait_profiles.png, contacts_all.npz, contacts_unique.npz
Usage: hic_remap_verify.py REF.fa PAIRS.tsv OUT [--bin_kb 500]
"""
import argparse
import os
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

JOINS = {"chr1_hap1": [263.0], "chr2_hap2": [215.0], "chr4_hap1": [85.5]}   # Mb
ARMS = {"chr1_hap1": "L1\u00b7P", "chr1_hap2": "L1\u00b7Q", "chr2_hap1": "L2\u00b7P", "chr2_hap2": "L2\u00b7Q"}
TESTS = [  # label, (left chrom, join), (right chrom, join), what the pollen said
    ("chr1_hap1: L1 | P   (same scaffold)", ("chr1_hap1", 263.0), ("chr1_hap1", 263.0), "pollen: not joined"),
    ("chr2_hap2: L2 | Q   (same scaffold)", ("chr2_hap2", 215.0), ("chr2_hap2", 215.0), "pollen: not joined"),
    ("L1 end | Q start   (across scaffolds)", ("chr1_hap1", 263.0), ("chr2_hap2", 215.0), "pollen: joined"),
    ("L2 end | P start   (across scaffolds)", ("chr2_hap2", 215.0), ("chr1_hap1", 263.0), "pollen: joined"),
    ("chr4_hap1: 0-85 Mb | rest   (chr3/chr4)", ("chr4_hap1", 85.5), ("chr4_hap1", 85.5), "pollen: joined"),
]
CONTROLS = [("chr1_hap1", 100), ("chr1_hap1", 180), ("chr1_hap1", 330), ("chr2_hap2", 100), ("chr2_hap2", 270),
            ("chr3_hap1", 150), ("chr4_hap1", 250), ("chr5_hap1", 60), ("chr6_hap1", 300)]
BACKGROUND = [(("chr1_hap1", 263.0), ("chr5_hap1", 100)), (("chr2_hap2", 215.0), ("chr6_hap1", 200)),
              (("chr3_hap1", 150), ("chr6_hap1", 300)), (("chr4_hap1", 250), ("chr5_hap1", 200))]
BAITS = [("L1 end  (chr1_hap1 257-262 Mb)", "chr1_hap1", 257, 262), ("P start  (chr1_hap1 264-269 Mb)", "chr1_hap1", 264, 269),
         ("L2 end  (chr2_hap2 209-214 Mb)", "chr2_hap2", 209, 214), ("Q start  (chr2_hap2 216-221 Mb)", "chr2_hap2", 216, 221)]
C12 = ["chr1_hap1", "chr1_hap2", "chr2_hap1", "chr2_hap2"]
BOX, GAP = 5.0, 1.0
NOTE = ("dashed: joins in the crossover reference. Real join: the red runs straight across. Misjoin: it pinches there, and "
        "each piece lights\nup against another scaffold instead. A faint copy of the diagonal between hap1 and hap2 = reads "
        "that fit either copy.")


class Genome:
    def __init__(self, lens, bin_bp):
        self.lens, self.bin, self.chroms = lens, bin_bp, list(lens)
        self.nb = {c: -(-lens[c] // bin_bp) for c in self.chroms}
        self.off, o = {}, 0
        for c in self.chroms:
            self.off[c] = o
            o += self.nb[c]
        self.n = o

    def rng(self, c, a_mb, b_mb):
        i0 = min(max(int(np.floor(a_mb * 1e6 / self.bin)), 0), self.nb[c])
        i1 = min(max(int(np.ceil(b_mb * 1e6 / self.bin)), 0), self.nb[c])
        return np.arange(self.off[c] + i0, self.off[c] + i1)


def count(path, G):
    Ma = np.zeros((G.n, G.n), dtype=np.float32)
    Mu = np.zeros((G.n, G.n), dtype=np.float32)
    st = defaultdict(lambda: np.zeros(6))          # library -> pairs, <2kb, 2kb-1Mb, >1Mb, other chrom, unique
    cols = ["c1", "p1", "q1", "c2", "p2", "q2", "lib"]
    for ch in pd.read_csv(path, sep="\t", header=None, usecols=[1, 2, 3, 4, 5, 6, 7], names=cols,
                          dtype={"c1": str, "p1": np.int64, "q1": np.int64, "c2": str, "p2": np.int64, "q2": np.int64,
                                 "lib": str}, chunksize=10_000_000):
        ch = ch[ch.c1.isin(G.off) & ch.c2.isin(G.off)]
        same = (ch.c1 == ch.c2).to_numpy()
        d = (ch.p1 - ch.p2).abs().to_numpy()
        uniq = ((ch.q1 >= 30) & (ch.q2 >= 30)).to_numpy()
        lib = ch.lib.to_numpy()
        for L in np.unique(lib):
            k = lib == L
            st[L] += [k.sum(), (k & same & (d < 2000)).sum(), (k & same & (d >= 2000) & (d < 1e6)).sum(),
                      (k & same & (d >= 1e6)).sum(), (k & ~same).sum(), (k & uniq & ~(same & (d < 2000))).sum()]
        keep = ~(same & (d < 2000))
        i = ch.c1.map(G.off).to_numpy(np.int64) + np.minimum(ch.p1.to_numpy() // G.bin, ch.c1.map(G.nb).to_numpy() - 1)
        j = ch.c2.map(G.off).to_numpy(np.int64) + np.minimum(ch.p2.to_numpy() // G.bin, ch.c2.map(G.nb).to_numpy() - 1)
        lo, hi = np.minimum(i, j), np.maximum(i, j)
        for M, sel in ((Ma, keep), (Mu, keep & uniq)):
            M += np.bincount(lo[sel] * G.n + hi[sel], minlength=G.n * G.n).reshape(G.n, G.n).astype(np.float32)
    Ma = Ma + Ma.T - np.diag(np.diag(Ma))
    Mu = Mu + Mu.T - np.diag(np.diag(Mu))
    return Ma, Mu, st


def box(M, G, left, right):
    (cl, jl), (cr, jr) = left, right
    I, J = G.rng(cl, jl - GAP - BOX, jl - GAP), G.rng(cr, jr + GAP, jr + GAP + BOX)
    return float(M[np.ix_(I, J)].sum()), len(I) * len(J)


def join_scores(M, G, name):
    P = G.off
    ctrl = [box(M, G, (c, x), (c, x)) for c, x in CONTROLS if c in P]
    bgv = [box(M, G, a, b) for a, b in BACKGROUND if a[0] in P and b[0] in P]
    exp = np.median([s / m for s, m in ctrl])
    bg = np.median([s / m for s, m in bgv])
    out = ["%s contacts: ordinary stretch of chromosome %.2f per bin pair (box totals %s); unrelated chromosomes %.3f"
           " = score %.3f (box totals %s)" % (name, exp, ", ".join("%d" % s for s, m in ctrl), bg, bg / max(exp, 1e-12),
                                              ", ".join("%d" % s for s, m in bgv))]
    for lab, a, b, pol in TESTS:
        if a[0] not in P or b[0] not in P:
            continue
        s, m = box(M, G, a, b)
        sc = s / m / max(exp, 1e-12)
        v = "joined" if sc >= 0.3 else ("not joined" if sc <= max(3 * bg / max(exp, 1e-12), 0.05) else "weak")
        if exp * m < 50:
            v = "too few contacts"
        out.append("  %-42s %9d contacts  score %6.3f   Hi-C: %-16s %s" % (lab, s, sc, v, pol))
    return out


def window_map(M, G, windows, agg, path, title):
    idx, gid, edges, g, labs = [], [], [0], 0, []
    for c, a, b in windows:
        I = G.rng(c, a, b)
        k = np.arange(len(I)) // agg
        idx.append(I)
        gid.append(g + k)
        g += int(k.max()) + 1 if len(k) else 0
        edges.append(g)
        labs.append("%s%s\n%g-%g Mb" % (c, " " + ARMS[c] if c in ARMS else "", a, b))
    idx, gid = np.concatenate(idx), np.concatenate(gid)
    P = np.zeros((len(idx), g), dtype=np.float32)
    P[np.arange(len(idx)), gid] = 1
    S = P.T @ M[np.ix_(idx, idx)] @ P
    L = np.log10(S + 1)
    pos = L[L > 0]
    fig, ax = plt.subplots(figsize=(9, 8.6))
    ax.imshow(L, cmap="Reds", vmin=0, vmax=np.percentile(pos, 99) if pos.size else 1, interpolation="nearest")
    for e in edges[1:-1]:
        ax.axhline(e - 0.5, color="#64748B", lw=0.7)
        ax.axvline(e - 0.5, color="#64748B", lw=0.7)
    mb = agg * G.bin / 1e6
    for k, (c, a, b) in enumerate(windows):
        for j in JOINS.get(c, []):
            if a <= j <= b:
                x = edges[k] + (j - a) / mb - 0.5
                ax.axvline(x, color="#0891B2", lw=1, ls="--")
                ax.axhline(x, color="#0891B2", lw=1, ls="--")
    mids = [(edges[k] + edges[k + 1]) / 2 - 0.5 for k in range(len(windows))]
    ax.set_xticks(mids)
    ax.set_xticklabels(labs, fontsize=8)
    ax.set_yticks(mids)
    ax.set_yticklabels(labs, fontsize=8)
    ax.set_title("%s  (%g Mb bins, %d contacts shown)" % (title, mb, S.sum() / 2), fontsize=10)
    ax.text(0, -0.09, NOTE, transform=ax.transAxes, fontsize=8, color="#334155", va="top")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def bait_profiles(M, G, path):
    cols = np.concatenate([G.rng(c, 0, G.lens[c] / 1e6) for c in C12])
    x = np.arange(len(cols))
    edges = np.cumsum([0] + [G.nb[c] for c in C12])
    fig, axs = plt.subplots(len(BAITS), 1, figsize=(12, 2.2 * len(BAITS)), sharex=True)
    for ax, (lab, c, a, b) in zip(axs, BAITS):
        y = M[np.ix_(G.rng(c, a, b), cols)].sum(axis=0)
        ax.bar(x, np.log10(y + 1), width=1.0, color="#534AB7")
        k = C12.index(c)
        ax.axvspan(edges[k] + a * 1e6 / G.bin, edges[k] + b * 1e6 / G.bin, color="#F59E0B", alpha=0.35, lw=0)
        for kk, cc in enumerate(C12):
            for j in JOINS.get(cc, []):
                ax.axvline(edges[kk] + j * 1e6 / G.bin, color="#0891B2", lw=1, ls="--")
        for e in edges[1:-1]:
            ax.axvline(e, color="#64748B", lw=1)
        ax.set_ylabel("log10\ncontacts", fontsize=8)
        ax.set_title("bait: %s  (orange)" % lab, fontsize=9, loc="left")
    axs[-1].set_xticks([(edges[k] + edges[k + 1]) / 2 for k in range(len(C12))])
    axs[-1].set_xticklabels(["%s  %s" % (c, ARMS[c]) for c in C12])
    fig.text(0.01, 0.0, "Where each arm end touches (unique reads). Besides its own neighbourhood, a piece lights up where it is "
             "physically joined: if L1 really continues into Q, the L1-end bait peaks at the start of Q on chr2_hap2.",
             fontsize=8, color="#334155", va="top")
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("pairs")
    ap.add_argument("out")
    ap.add_argument("--bin_kb", type=int, default=500)
    A = ap.parse_args()
    os.makedirs(A.out, exist_ok=True)
    lens = {}
    for l in open(A.ref + ".fai"):
        f = l.split("\t")
        lens[f[0]] = int(f[1])
    G = Genome(lens, A.bin_kb * 1000)
    Ma, Mu, st = count(A.pairs, G)
    rows = ["HI-C LIBRARIES  (read pairs with both mates placed, MAPQ >= 1)",
            "  %-44s %12s %9s %9s %9s %9s %12s" % ("library", "pairs", "<2 kb", "2kb-1Mb", ">1 Mb", "other chr", "unique, usable")]
    for L, v in sorted(st.items()) + [("ALL", sum(st.values()))]:
        n = max(v[0], 1)
        rows.append("  %-44s %12d %8.1f%% %8.1f%% %8.1f%% %8.1f%% %12d" % (
            L[:44], v[0], 100 * v[1] / n, 100 * v[2] / n, 100 * v[3] / n, 100 * v[4] / n, v[5]))
    rows.append("  (<2 kb = not a ligation; usable Hi-C = the other three columns)")
    open(os.path.join(A.out, "library_stats.txt"), "w").write("\n".join(rows) + "\n")
    print("\n".join(rows))
    js = ["HI-C JOIN TEST  (%g kb bins; %g Mb boxes either side of a join, %g Mb left out each side)" % (
        G.bin / 1e3, BOX, GAP)] + join_scores(Ma, G, "ALL (MAPQ >= 1)") + [""] + join_scores(Mu, G, "UNIQUE (MAPQ >= 30)")
    open(os.path.join(A.out, "join_scores.txt"), "w").write("\n".join(js) + "\n")
    print("\n".join(js))
    np.savez_compressed(os.path.join(A.out, "contacts_all.npz"), M=Ma, chroms=np.array(G.chroms), bin=G.bin)
    np.savez_compressed(os.path.join(A.out, "contacts_unique.npz"), M=Mu, chroms=np.array(G.chroms), bin=G.bin)
    agg = max(1, int(round(2e6 / G.bin)))
    whole = lambda c: (c, 0, int(lens[c] / 1e6) + 1)
    window_map(Ma, G, [whole(c) for c in C12 if c in lens], agg, os.path.join(A.out, "hicmap_chr1_chr2.png"),
               "chr1 and chr2, both haplotypes, all reads")
    window_map(Ma, G, [("chr1_hap1", 223, 303), ("chr2_hap2", 175, 255)], 1, os.path.join(A.out, "hicmap_chr1_chr2_zoom.png"),
               "around the two joins, all reads")
    window_map(Mu, G, [("chr1_hap1", 223, 303), ("chr2_hap2", 175, 255)], 1,
               os.path.join(A.out, "hicmap_chr1_chr2_zoom_unique.png"), "around the two joins, unique reads only")
    window_map(Ma, G, [("chr1_hap2", 222, 302), ("chr2_hap1", 171, 251)], 1,
               os.path.join(A.out, "hicmap_chr1_chr2_zoom_control.png"), "control: the same arms joined the other way, all reads")
    for a, b in (("chr3", "chr4"), ("chr5", "chr6")):
        cs = [c for c in ("%s_hap1" % a, "%s_hap2" % a, "%s_hap1" % b, "%s_hap2" % b) if c in lens]
        window_map(Ma, G, [whole(c) for c in cs], agg, os.path.join(A.out, "hicmap_%s_%s.png" % (a, b)),
                   "%s and %s, both haplotypes, all reads" % (a, b))
    bait_profiles(Mu, G, os.path.join(A.out, "bait_profiles.png"))
    print("wrote:\n  " + "\n  ".join(sorted(os.path.join(os.path.abspath(A.out), f) for f in os.listdir(A.out)
                                             if f.endswith((".png", ".txt")))))


if __name__ == "__main__":
    main()
