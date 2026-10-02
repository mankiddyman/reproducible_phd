#!/usr/bin/env python3
"""hic_verify.py -- do the Hi-C contacts support each join in the D. paradoxa assembly?

Hi-C read pairs from the same plant are mapped to the final dual assembly (both
haplotypes, 12 chromosomes), each mate on its own, keeping pairs where both mates
have MAPQ >= 30. Input, one line per pair: name, chrom1, pos1, chrom2, pos2.

Pieces joined in the nucleus contact each other like any stretch of one
chromosome; pieces on different chromosomes contact only at background level.
Join score = mean contacts between the 5 Mb either side of a join (1 Mb left out
on each side), divided by the same-shaped box inside ordinary chromosome arms
(median of controls). ~1: joined in the nucleus. Down at the unrelated-
chromosome level: not joined.

Writes to OUT:
  hic_chr1_chr2.png, hic_chr3_chr4.png, hic_chr5_chr6.png  both haplotypes, 2 Mb bins
  hic_join_scores.png / .txt   the join test as bars, next to what the pollen said
  hic_vs_pollen.png            Hi-C and pollen on the crossover reference, same windows
  agp_at_joins.txt             what the assembly put at each join: contigs, gaps, RagTag units
  contacts.npz                 the contact matrix

Usage: hic_verify.py REF.fa PAIRS.tsv OUT [--bin_kb 500] [--agp FINAL.agp]
                     [--ragtag RAGTAG.agp] [--pollen LINKAGE_DIR]
"""
import argparse
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Joins under test, in Mb. chr1/chr2: P starts at 265 on chr1_hap1 (gap at 262.93),
# Q starts at 215 on chr2_hap2. chr4 0-85.5 Mb is the chr3/chr4 translocated piece.
TESTS = [  # label, (left chrom, join), (right chrom, join), what the pollen said
    ("chr1_hap1: L1 | P   (same scaffold)", ("chr1_hap1", 264.0), ("chr1_hap1", 264.0), "pollen: not joined"),
    ("chr2_hap2: L2 | Q   (same scaffold)", ("chr2_hap2", 216.0), ("chr2_hap2", 216.0), "pollen: not joined"),
    ("L1 end | Q start   (across scaffolds)", ("chr1_hap1", 264.0), ("chr2_hap2", 216.0), "pollen: joined"),
    ("L2 end | P start   (across scaffolds)", ("chr2_hap2", 216.0), ("chr1_hap1", 264.0), "pollen: joined"),
    ("chr4_hap1: 0-85 Mb | rest   (chr3/chr4)", ("chr4_hap1", 85.5), ("chr4_hap1", 85.5), "pollen: joined"),
]
CONTROLS = [("chr1_hap1", 100), ("chr1_hap1", 180), ("chr1_hap1", 330), ("chr2_hap2", 100),
            ("chr2_hap2", 270), ("chr3_hap1", 150), ("chr4_hap1", 250), ("chr5_hap1", 60),
            ("chr6_hap1", 300)]
BACKGROUND = [(("chr1_hap1", 264.0), ("chr5_hap1", 100)), (("chr2_hap2", 216.0), ("chr6_hap1", 200)),
              (("chr3_hap1", 150), ("chr6_hap1", 300)), (("chr4_hap1", 250), ("chr5_hap1", 200))]
MARKS = {"chr1_hap1": [264.0], "chr2_hap2": [216.0], "chr4_hap1": [85.5]}
ARMS = {"chr1_hap1": "L1\u00b7P", "chr1_hap2": "L1\u00b7Q", "chr2_hap1": "L2\u00b7P", "chr2_hap2": "L2\u00b7Q"}
PAIRS = [("chr1", "chr2"), ("chr3", "chr4"), ("chr5", "chr6")]
AGP_WINDOWS = [("chr1_hap1", 255, 272), ("chr2_hap2", 203, 228), ("chr1_hap2", 250, 280),
               ("chr2_hap1", 200, 225), ("chr4_hap1", 80, 92)]
BOX, GAP = 5.0, 1.0      # Mb
AGG = 4                  # map bins: 4 x 500 kb = 2 Mb


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
        """bin indices covering [a, b) Mb of chromosome c, clipped to its length"""
        i0 = int(np.floor(a_mb * 1e6 / self.bin))
        i1 = int(np.ceil(b_mb * 1e6 / self.bin))
        i0, i1 = min(max(i0, 0), self.nb[c]), min(max(i1, 0), self.nb[c])
        return np.arange(self.off[c] + i0, self.off[c] + i1)


def read_fai(fa):
    lens = {}
    for l in open(fa + ".fai"):
        f = l.split("\t")
        lens[f[0]] = int(f[1])
    return lens


def count(path, G):
    M = np.zeros((G.n, G.n), dtype=np.float32)
    st = dict(used=0, close=0, other=0)
    for ch in pd.read_csv(path, sep="\t", header=None, usecols=[1, 2, 3, 4], names=["c1", "p1", "c2", "p2"],
                          dtype={"c1": str, "p1": np.int64, "c2": str, "p2": np.int64}, chunksize=10_000_000):
        keep = ch.c1.isin(G.off) & ch.c2.isin(G.off)
        st["other"] += int((~keep).sum())
        ch = ch[keep]
        near = (ch.c1 == ch.c2) & ((ch.p1 - ch.p2).abs() < 2000)
        st["close"] += int(near.sum())
        ch = ch[~near]
        i = ch.c1.map(G.off).to_numpy(np.int64) + np.minimum(ch.p1.to_numpy() // G.bin, ch.c1.map(G.nb).to_numpy() - 1)
        j = ch.c2.map(G.off).to_numpy(np.int64) + np.minimum(ch.p2.to_numpy() // G.bin, ch.c2.map(G.nb).to_numpy() - 1)
        lo, hi = np.minimum(i, j), np.maximum(i, j)
        M += np.bincount(lo * G.n + hi, minlength=G.n * G.n).reshape(G.n, G.n).astype(np.float32)
        st["used"] += len(ch)
    M = M + M.T - np.diag(np.diag(M))
    return M, st


def box(M, G, left, right):
    (cl, jl), (cr, jr) = left, right
    I = G.rng(cl, jl - GAP - BOX, jl - GAP)
    J = G.rng(cr, jr + GAP, jr + GAP + BOX)
    return float(M[np.ix_(I, J)].mean()) if len(I) and len(J) else np.nan


def join_scores(M, G, out):
    ctrl = [box(M, G, (c, x), (c, x)) for c, x in CONTROLS if c in G.off]
    bgv = [box(M, G, a, b) for a, b in BACKGROUND if a[0] in G.off and b[0] in G.off]
    exp, bg = np.nanmedian(ctrl), np.nanmedian(bgv)
    L = ["HI-C JOIN TEST  (mean contacts per %g kb x %g kb bin pair; %g Mb boxes either side of the join, %g Mb left out each side)" % (
             G.bin / 1e3, G.bin / 1e3, BOX, GAP),
         "  score = contacts / an ordinary stretch of chromosome; ~1 = joined in the nucleus, ~background = not joined",
         "  ordinary stretch of chromosome (median of %d controls): %.2f   [%s]" % (
             len(ctrl), exp, ", ".join("%.2f" % v for v in ctrl)),
         "  unrelated chromosomes (median of %d): %.2f = score %.3f" % (len(bgv), bg, bg / exp)]
    rows = []
    for lab, a, b, pol in TESTS:
        if a[0] not in G.off or b[0] not in G.off:
            continue
        v = box(M, G, a, b)
        s = v / exp
        hic = "joined" if s >= 0.3 else ("not joined" if s <= max(3 * bg / exp, 0.05) else "weak")
        rows.append((lab, s, hic, pol))
        L.append("  %-42s %8.2f  score %.3f   Hi-C: %-10s  %s" % (lab, v, s, hic, pol))
    open(os.path.join(out, "hic_join_scores.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))
    fig, ax = plt.subplots(figsize=(9, 4.2))
    labs = [r[0] for r in rows] + ["unrelated chromosomes"]
    vals = [r[1] for r in rows] + [bg / exp]
    y = np.arange(len(labs))[::-1]
    ax.barh(y, vals, color=["#534AB7" if v >= 0.3 else "#CECBF6" for v in vals], height=0.6)
    ax.axvline(1, color="#334155", lw=1, ls="--")
    ax.text(1.02, y[0] + 0.45, "ordinary stretch\nof chromosome", fontsize=8, color="#334155", va="top")
    for yy, (lab, s, hic, pol) in zip(y, rows):
        ax.text(max(s + 0.03, 1.05) if s > 0.6 else s + 0.03, yy, "Hi-C: %s   (%s)" % (hic, pol),
                va="center", fontsize=8, color="#334155")
    ax.text(bg / exp + 0.03, y[-1], "background", va="center", fontsize=8, color="#334155")
    ax.set_yticks(y)
    ax.set_yticklabels(labs, fontsize=9)
    ax.set_xlim(0, max(1.6, np.nanmax(vals) * 1.25 + 0.6))
    ax.set_xlabel("Hi-C contacts across the join, relative to an ordinary stretch of chromosome")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "hic_join_scores.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def onehot(groups, n_groups):
    P = np.zeros((len(groups), n_groups), dtype=np.float32)
    P[np.arange(len(groups)), groups] = 1
    return P


def pair_map(M, G, cs, path, title):
    idx, gid, edges, g = [], [], [0], 0
    for c in cs:
        k = np.arange(G.nb[c]) // AGG
        idx.append(np.arange(G.off[c], G.off[c] + G.nb[c]))
        gid.append(g + k)
        g += int(k.max()) + 1
        edges.append(g)
    idx, gid = np.concatenate(idx), np.concatenate(gid)
    P = onehot(gid, g)
    S = P.T @ M[np.ix_(idx, idx)] @ P
    L = np.log10(S + 1)
    vmax = np.percentile(L[L > 0], 99.5) if (L > 0).any() else 1
    fig, ax = plt.subplots(figsize=(9, 8.6))
    ax.imshow(L, cmap="Reds", vmin=0, vmax=vmax, interpolation="nearest")
    for e in edges[1:-1]:
        ax.axhline(e - 0.5, color="#64748B", lw=0.6)
        ax.axvline(e - 0.5, color="#64748B", lw=0.6)
    mb = AGG * G.bin / 1e6
    for k, c in enumerate(cs):
        for j in MARKS.get(c, []):
            x = edges[k] + j / mb - 0.5
            ax.axvline(x, color="#0891B2", lw=1, ls="--")
            ax.axhline(x, color="#0891B2", lw=1, ls="--")
    mids = [(edges[k] + edges[k + 1]) / 2 - 0.5 for k in range(len(cs))]
    labs = [c + ("\n" + ARMS[c] if c in ARMS else "") for c in cs]
    ax.set_xticks(mids)
    ax.set_xticklabels(labs, fontsize=9)
    ax.set_yticks(mids)
    ax.set_yticklabels(labs, fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.text(0, -0.08, "dashed: a join in the current crossover reference. Real join: the red continues across it. "
            "Misjoin: it stops there,\nand the piece lights up against another scaffold instead. "
            "Faint copy of the diagonal between hap1 and hap2 = reads that fit either copy.",
            transform=ax.transAxes, fontsize=8, color="#334155", va="top")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def hic_vs_pollen(M, G, pdir, path):
    win = pd.read_csv(os.path.join(pdir, "windows.tsv"), sep="\t")
    R = np.load(os.path.join(pdir, "r_matrix.npy"))
    W = 10.0
    s = os.path.join(pdir, "linkage_summary.txt")
    if os.path.exists(s):
        m = re.search(r"([\d.]+) Mb windows", open(s).read())
        if m:
            W = float(m.group(1))
    chroms = list(dict.fromkeys(win.chrom))
    rows, cols, k0, edges = [], [], 0, [0]
    for c in chroms:
        nw = int((win.chrom == c).sum())
        b = np.arange(G.nb[c])
        rows.append(G.off[c] + b)
        cols.append(k0 + np.minimum((b * G.bin) // int(W * 1e6), nw - 1).astype(int))
        k0 += nw
        edges.append(k0)
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    P = onehot(cols, k0)
    H = P.T @ M[np.ix_(rows, rows)] @ P
    LH = np.log10(H + 1)
    fig, axs = plt.subplots(1, 2, figsize=(15, 7.4))
    axs[0].imshow(LH, cmap="magma_r", vmin=np.percentile(LH, 5), vmax=np.percentile(LH, 99.5), interpolation="nearest")
    axs[0].set_title("Hi-C: physical contact (dark = touching)", fontsize=11)
    axs[1].imshow(R, cmap="magma", vmin=0, vmax=0.5, interpolation="nearest")
    axs[1].set_title("Pollen: inherited together (dark = linked, r near 0)", fontsize=11)
    mids = [(edges[k] + edges[k + 1]) / 2 - 0.5 for k in range(len(chroms))]
    for ax in axs:
        for e in edges[1:-1]:
            ax.axhline(e - 0.5, color="white", lw=0.8)
            ax.axvline(e - 0.5, color="white", lw=0.8)
        ax.set_xticks(mids)
        ax.set_xticklabels(chroms, rotation=90, fontsize=8)
        ax.set_yticks(mids)
        ax.set_yticklabels(chroms, fontsize=8)
    fig.text(0.01, 0.0, "Same reference, same %g Mb windows. A misjoin shows in both maps: the piece really sits elsewhere.\n"
             "A real translocation heterozygote shows only in pollen: its pieces travel together through meiosis, "
             "but sit on separate molecules within each haplotype." % W, fontsize=9, color="#334155", va="top")
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def fetch(fa, idx, c, s, e):
    L, off, lb, lw = idx[c]
    s, e = max(0, s), min(L, e)
    b0 = off + (s // lb) * lw + s % lb
    b1 = off + (e // lb) * lw + e % lb
    with open(fa, "rb") as fh:
        fh.seek(b0)
        return fh.read(b1 - b0).replace(b"\n", b"").replace(b"\r", b"")


def read_agp(path):
    rows = {}
    for l in open(path):
        if l.startswith("#") or not l.strip():
            continue
        f = l.rstrip("\n").split("\t")
        rows.setdefault(f[0], []).append(f)
    return rows


def agp_report(agp, ragtag, ref, lens, out_path):
    idx = {}
    for l in open(ref + ".fai"):
        f = l.split("\t")
        idx[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    fin = read_agp(agp)
    rag = read_agp(ragtag) if ragtag and os.path.exists(ragtag) else {}
    rag_of = {}
    for o, fs in rag.items():
        for f in fs:
            if f[4] not in ("N", "U"):
                rag_of[f[5]] = o
    by_len = {}
    for c, L in lens.items():
        by_len.setdefault(L, []).append(c)
    name = {}
    for s, fs in fin.items():
        L = max(int(f[2]) for f in fs)
        if len(by_len.get(L, [])) == 1:
            name[by_len[L][0]] = s
    out = ["WHAT THE ASSEMBLY PUT AROUND EACH JOIN",
           "  final AGP: %s" % agp, "  RagTag AGP: %s" % (ragtag or "-"),
           "  contig = made by the assembler from HiFi reads; gap between HapHiC units = joined by Hi-C",
           "  scaffolding or curation; gap inside a RagTag unit = joined by RagTag from a reference", ""]
    for c, a, b in AGP_WINDOWS:
        if c not in name:
            out.append("%s: no AGP object of the same length" % c)
            continue
        s, L, fs = name[c], lens[c], fin[name[c]]
        gaps = [(int(f[1]), int(f[5])) for f in fs if f[4] in ("N", "U")][:25]
        same = flip = 0
        for g0, gl in gaps:
            k = min(gl, 10)
            same += fetch(ref, idx, c, g0 - 1, g0 - 1 + k).upper() == b"N" * k
            q = L - (g0 - 1) - gl
            flip += fetch(ref, idx, c, q, q + k).upper() == b"N" * k
        rev = flip > same
        out.append("%s = %s, %s orientation (gaps matched: same %d, reversed %d)  window %g-%g Mb" % (
            c, s, "reversed" if rev else "same", same, flip, a, b))
        lines = []
        for f in fs:
            ob, oe = int(f[1]), int(f[2])
            cb, ce = (L - oe + 1, L - ob + 1) if rev else (ob, oe)
            if ce < a * 1e6 or cb > b * 1e6:
                continue
            if f[4] in ("N", "U"):
                lines.append((cb, "  %9.3f-%9.3f  gap %s bp  %s" % (cb / 1e6, ce / 1e6, f[5], " ".join(f[6:9]))))
                continue
            comp, xb, xe, o = f[5], int(f[6]), int(f[7]), f[8]
            if rev:
                o = {"+": "-", "-": "+"}.get(o, o)
            where = ("inside RagTag unit %s" % rag_of[comp]) if comp in rag_of else (
                "RagTag unit itself" if comp in rag else "not in the RagTag AGP")
            lines.append((cb, "  %9.3f-%9.3f  %-30s %s  %s" % (cb / 1e6, ce / 1e6, comp, o, where)))
            if comp in rag:     # show RagTag's own contigs and gaps inside this unit, in chromosome coordinates
                lo, hi = max(cb, a * 1e6), min(ce, b * 1e6)
                for g in rag[comp]:
                    rb, re_ = int(g[1]), int(g[2])
                    if o == "+":
                        gb, ge = cb + (rb - xb), cb + (re_ - xb)
                    else:
                        gb, ge = cb + (xe - re_), cb + (xe - rb)
                    if ge < lo or gb > hi:
                        continue
                    go = g[8] if o == "+" else {"+": "-", "-": "+"}.get(g[8], g[8])
                    what = ("RagTag gap %s bp %s" % (g[5], " ".join(g[6:9]))) if g[4] in ("N", "U") else (
                        "contig %s %s" % (g[5], go))
                    lines.append((gb + 0.5, "      %9.3f-%9.3f    %s" % (gb / 1e6, ge / 1e6, what)))
        out += [t for _, t in sorted(lines)] + [""]
    open(out_path, "w").write("\n".join(out) + "\n")
    print("\n".join(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("pairs")
    ap.add_argument("out")
    ap.add_argument("--bin_kb", type=int, default=500)
    ap.add_argument("--agp")
    ap.add_argument("--ragtag")
    ap.add_argument("--pollen")
    A = ap.parse_args()
    os.makedirs(A.out, exist_ok=True)
    lens = read_fai(A.ref)
    G = Genome(lens, A.bin_kb * 1000)
    M, st = count(A.pairs, G)
    print("pairs used %d | mates < 2 kb apart, dropped %d | off the 12 chromosomes %d" % (
        st["used"], st["close"], st["other"]))
    np.savez_compressed(os.path.join(A.out, "contacts.npz"), M=M, chroms=np.array(G.chroms),
                        nb=np.array([G.nb[c] for c in G.chroms]), bin=G.bin)
    join_scores(M, G, A.out)
    for a, b in PAIRS:
        cs = [c for c in ("%s_hap1" % a, "%s_hap2" % a, "%s_hap1" % b, "%s_hap2" % b) if c in G.off]
        pair_map(M, G, cs, os.path.join(A.out, "hic_%s_%s.png" % (a, b)),
                 "D. paradoxa Hi-C, %s and %s, both haplotypes (%g Mb bins, %d pairs)" % (
                     a, b, AGG * G.bin / 1e6, st["used"]))
    if A.pollen and os.path.exists(os.path.join(A.pollen, "r_matrix.npy")):
        hic_vs_pollen(M, G, A.pollen, os.path.join(A.out, "hic_vs_pollen.png"))
    if A.agp and os.path.exists(A.agp):
        agp_report(A.agp, A.ragtag, A.ref, lens, os.path.join(A.out, "agp_at_joins.txt"))
    print("wrote:\n  " + "\n  ".join(sorted(os.path.join(os.path.abspath(A.out), f) for f in os.listdir(A.out))))


if __name__ == "__main__":
    main()
