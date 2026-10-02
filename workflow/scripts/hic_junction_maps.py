#!/usr/bin/env python3
"""hic_junction_maps.py -- the full-depth Hi-C contact matrix around the D. paradoxa joins.

Reads the Juicebox map HapHiC built in pass 2 from all the Hi-C reads
(04.build/out_JBAT.hic; its coordinates are HapHiC's units laid end to end in the
order of out_JBAT.assembly, divided by an integer scale factor) and redraws it in
the coordinates of the final chromosomes, using the final AGP (final chromosome
-> unit, position, orientation). No re-mapping: these are the contacts you curated on.

Usage: hic_junction_maps.py OUT_JBAT.hic OUT_JBAT.assembly FINAL.agp CHR.fa.fai OUTDIR [HIC_CHROM]
Writes to OUTDIR: hicmap_chr1_chr2.png, hicmap_chr1_chr2_zoom.png, hicmap_chr1_chr2_zoom_control.png,
                  hicmap_chr3_chr4.png, hicmap_chr5_chr6.png, hic_join_scores_fulldepth.txt
"""
import os
import sys
from collections import defaultdict

import numpy as np
import hicstraw
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

JOINS = {"chr1_hap1": [263.0], "chr2_hap2": [215.0], "chr4_hap1": [85.5]}   # Mb, final coordinates
ARMS = {"chr1_hap1": "L1\u00b7P", "chr1_hap2": "L1\u00b7Q", "chr2_hap1": "L2\u00b7P", "chr2_hap2": "L2\u00b7Q"}
TESTS = [  # label, (left chrom, join Mb), (right chrom, join Mb), what the pollen said
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
BOX, GAP = 5.0, 1.0      # Mb: 5 Mb either side of a join, 1 Mb left out on each side
NOTE = ("dashed: joins in the crossover reference. Real join: the red runs straight across. Misjoin: it pinches there, "
        "and each piece lights\nup against another scaffold instead. A faint copy of the diagonal between hap1 and hap2 "
        "= reads that fit either copy.")
G = {}                   # filled by setup(): hic, chrom, K, RES, unit_at, pieces, chr_len


def setup(hic_path, asm_path, agp_path, fai_path, hic_chrom="assembly"):
    frag, order = {}, []
    for l in open(asm_path):
        l = l.strip()
        if not l:
            continue
        if l.startswith(">"):
            f = l[1:].split()
            frag[int(f[1])] = (f[0], int(f[2]))
        else:
            order.append([int(x) for x in l.split()])
    unit_at, cum = {}, 0                  # unit -> (start in assembly coordinates, length, orientation)
    for line in order:
        for s in line:
            name, ln = frag[abs(s)]
            unit_at[name] = (cum, ln, "+" if s > 0 else "-")
            cum += ln
    hic = hicstraw.HiCFile(hic_path)
    hic_len = {c.name: c.length for c in hic.getChromosomes()}
    if hic_chrom not in hic_len:
        sys.exit("no chromosome %s in %s: %s" % (hic_chrom, hic_path, list(hic_len)[:5]))
    K = max(1, int(round(cum / hic_len[hic_chrom])))
    print("JBAT: %d units, %.3f Gb end to end; .hic '%s' = %.3f Gb -> scale factor %d (check %.4f, should be ~1)" % (
        len(unit_at), cum / 1e9, hic_chrom, hic_len[hic_chrom] / 1e9, K, cum / K / hic_len[hic_chrom]))
    chrom_by_len = {}
    for l in open(fai_path):
        f = l.split("\t")
        chrom_by_len[int(f[1])] = f[0]
    rows = defaultdict(list)
    for l in open(agp_path):
        if l.startswith("#") or not l.strip():
            continue
        f = l.rstrip("\n").split("\t")
        rows[f[0]].append(f)
    pieces = {}
    for s, fs in rows.items():
        c = chrom_by_len.get(max(int(f[2]) for f in fs))
        if c:
            pieces[c] = [(int(f[1]), int(f[2]), f[5], int(f[6]), int(f[7]), f[8]) for f in fs if f[4] not in ("N", "U")]
    missing = sorted({p[2] for ps in pieces.values() for p in ps if p[2] not in unit_at})
    print("final chromosomes found: %d; final-AGP units missing from out_JBAT.assembly: %s" % (
        len(pieces), ", ".join(missing[:8]) if missing else "none"))
    G.update(hic=hic, chrom=hic_chrom, K=K, RES=sorted(hic.getResolutions()), unit_at=unit_at, pieces=pieces,
             chr_len={v: k for k, v in chrom_by_len.items()})


def to_asm(c, p):
    """final chromosome c, position p (1-based) -> JBAT assembly coordinate (0-based), or None"""
    for ob, oe, comp, cb, ce, o in G["pieces"][c]:
        if ob <= p <= oe and comp in G["unit_at"]:
            q = cb + (p - ob) if o != "-" else ce - (p - ob)
            s, ln, oj = G["unit_at"][comp]
            return s + (q - 1) if oj == "+" else s + (ln - q)
    return None


def track(windows, r):
    """windows [(chrom, start Mb, end Mb)] -> segments of contiguous .hic bins at resolution r (.hic units)"""
    segs, n, ticks = [], 0, []
    for c, a, b in windows:
        x0, x1 = int(a * 1e6) + 1, int(b * 1e6)
        t0 = n
        for ob, oe, comp, cb, ce, o in G["pieces"][c]:
            lo, hi = max(ob, x0), min(oe, x1)
            if lo > hi or comp not in G["unit_at"]:
                continue
            b_lo, b_hi = to_asm(c, lo) // G["K"] // r, to_asm(c, hi) // G["K"] // r
            m = abs(b_hi - b_lo) + 1
            segs.append(dict(i0=n, n=m, b_lo=b_lo, b_hi=b_hi, step=1 if b_hi >= b_lo else -1,
                             chrom=c, mb0=lo / 1e6, mb1=hi / 1e6))
            n += m
        ticks.append((c, a, b, t0, n))
    return segs, n, ticks


def matrix(seg1, n1, seg2, n2, r):
    mzd = G["hic"].getMatrixZoomData(G["chrom"], G["chrom"], "observed", "NONE", "BP", r)
    M = np.zeros((n1, n2))
    for s in seg1:
        for t in seg2:
            x0, x1 = sorted((s["b_lo"], s["b_hi"]))
            y0, y1 = sorted((t["b_lo"], t["b_hi"]))
            B = np.asarray(mzd.getRecordsAsMatrix(x0 * r, x1 * r, y0 * r, y1 * r), dtype=float)
            if B.shape != (x1 - x0 + 1, y1 - y0 + 1):
                B2 = np.zeros((x1 - x0 + 1, y1 - y0 + 1))
                B2[:min(B.shape[0], B2.shape[0]), :min(B.shape[1], B2.shape[1])] = \
                    B[:B2.shape[0], :B2.shape[1]]
                B = B2
            if s["step"] < 0:
                B = B[::-1, :]
            if t["step"] < 0:
                B = B[:, ::-1]
            M[s["i0"]:s["i0"] + s["n"], t["i0"]:t["i0"] + t["n"]] += B
    return M


def pick(real_bp):
    return min(G["RES"], key=lambda r: abs(np.log(r * G["K"] / real_bp)))


def bin_of(segs, c, mb):
    for s in segs:
        if s["chrom"] == c and min(s["mb0"], s["mb1"]) <= mb <= max(s["mb0"], s["mb1"]):
            return s["i0"] + (mb - s["mb0"]) / max(s["mb1"] - s["mb0"], 1e-9) * (s["n"] - 1)
    return None


def draw(windows, real_bp, path, title, note=NOTE):
    r = pick(real_bp)
    segs, n, ticks = track(windows, r)
    M = matrix(segs, n, segs, n, r)
    I, J = np.indices(M.shape)
    d = np.abs(I - J)
    near = M[d == 1].mean() if (d == 1).any() else np.nan
    far = M[d >= 20].mean() if (d >= 20).any() else np.nan
    print("%-34s %4d bins of %4.0f kb, %9d contacts; next to diagonal %.1f vs 20+ bins away %.3f per bin pair" % (
        os.path.basename(path), n, r * G["K"] / 1e3, M.sum() / 2, near, far))
    fig, ax = plt.subplots(figsize=(9, 8.6))
    L = np.log10(M + 1)
    pos = L[L > 0]
    ax.imshow(L, cmap="Reds", vmin=0, vmax=np.percentile(pos, 99) if pos.size else 1, interpolation="nearest")
    for c, a, b, t0, t1 in ticks[1:]:
        ax.axhline(t0 - 0.5, color="#64748B", lw=0.7)
        ax.axvline(t0 - 0.5, color="#64748B", lw=0.7)
    for c, js in JOINS.items():
        for j in js:
            x = bin_of(segs, c, j)
            if x is not None:
                ax.axvline(x, color="#0891B2", lw=1, ls="--")
                ax.axhline(x, color="#0891B2", lw=1, ls="--")
    mids = [(t0 + t1) / 2 - 0.5 for c, a, b, t0, t1 in ticks]
    labs = ["%s%s\n%g-%g Mb" % (c, " " + ARMS[c] if c in ARMS else "", a, b) for c, a, b, t0, t1 in ticks]
    ax.set_xticks(mids)
    ax.set_xticklabels(labs, fontsize=8)
    ax.set_yticks(mids)
    ax.set_yticklabels(labs, fontsize=8)
    ax.set_title("%s  (%.0f kb bins, full Hi-C from HapHiC pass 2)" % (title, r * G["K"] / 1e3), fontsize=10)
    ax.text(0, -0.09, note, transform=ax.transAxes, fontsize=8, color="#334155", va="top")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return M


def join_test(path):
    r = pick(250e3)

    def box(left, right):
        (cl, jl), (cr, jr) = left, right
        s1, n1, _ = track([(cl, jl - GAP - BOX, jl - GAP)], r)
        s2, n2, _ = track([(cr, jr + GAP, jr + GAP + BOX)], r)
        M = matrix(s1, n1, s2, n2, r)
        return M.sum(), M.size

    P = G["pieces"]
    ctrl = [box((c, x), (c, x)) for c, x in CONTROLS if c in P]
    bgv = [box(a, b) for a, b in BACKGROUND if a[0] in P and b[0] in P]
    exp = np.median([s / m for s, m in ctrl])
    bg = np.median([s / m for s, m in bgv])
    out = ["HI-C JOIN TEST, FULL DEPTH  (HapHiC pass-2 .hic; %.0f kb bins; %g Mb boxes either side of a join, %g Mb left out)" % (
               r * G["K"] / 1e3, BOX, GAP),
           "  ordinary stretch of chromosome: %.2f contacts per bin pair (median of %d; box totals %s)" % (
               exp, len(ctrl), ", ".join("%d" % s for s, m in ctrl)),
           "  unrelated chromosomes: %.3f per bin pair = score %.3f (box totals %s)" % (
               bg, bg / exp, ", ".join("%d" % s for s, m in bgv)), ""]
    for lab, a, b, pol in TESTS:
        if a[0] not in P or b[0] not in P:
            continue
        s, m = box(a, b)
        sc = s / m / exp
        v = "joined" if sc >= 0.3 else ("not joined" if sc <= max(3 * bg / exp, 0.05) else "weak")
        if exp * m < 50:
            v = "too few contacts"
        out.append("  %-42s %8d contacts  score %.3f   Hi-C: %-16s %s" % (lab, s, sc, v, pol))
    open(path, "w").write("\n".join(out) + "\n")
    print("\n".join(out))


def main():
    hic_path, asm_path, agp_path, fai_path, out = sys.argv[1:6]
    setup(hic_path, asm_path, agp_path, fai_path, sys.argv[6] if len(sys.argv) > 6 else "assembly")
    os.makedirs(out, exist_ok=True)
    P, CL = G["pieces"], G["chr_len"]
    whole = lambda c: (c, 0, int(CL[c] / 1e6) + 1)
    draw([whole(c) for c in ("chr1_hap1", "chr1_hap2", "chr2_hap1", "chr2_hap2") if c in P], 2e6,
         os.path.join(out, "hicmap_chr1_chr2.png"), "chr1 and chr2, both haplotypes")
    draw([("chr1_hap1", 223, 303), ("chr2_hap2", 175, 255)], 250e3,
         os.path.join(out, "hicmap_chr1_chr2_zoom.png"), "around the two joins in question")
    draw([("chr1_hap2", 222, 302), ("chr2_hap1", 171, 251)], 250e3,
         os.path.join(out, "hicmap_chr1_chr2_zoom_control.png"), "the other pair: same arms, joined the other way round")
    draw([whole(c) for c in ("chr3_hap1", "chr3_hap2", "chr4_hap1", "chr4_hap2") if c in P], 2e6,
         os.path.join(out, "hicmap_chr3_chr4.png"), "chr3 and chr4, both haplotypes")
    draw([whole(c) for c in ("chr5_hap1", "chr5_hap2", "chr6_hap1", "chr6_hap2") if c in P], 2e6,
         os.path.join(out, "hicmap_chr5_chr6.png"), "chr5 and chr6, both haplotypes")
    join_test(os.path.join(out, "hic_join_scores_fulldepth.txt"))
    print("wrote:\n  " + "\n  ".join(sorted(os.path.join(os.path.abspath(out), f) for f in os.listdir(out)
                                             if f.startswith(("hicmap", "hic_join_scores_fulldepth")))))


if __name__ == "__main__":
    main()
