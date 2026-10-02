#!/usr/bin/env python3
"""join_provenance.py -- which step of the assembly made each suspect join?

The joins sit on RagTag 'align_genus' gaps: pass-2 RagTag laid the pass-2 contigs
onto its reference (the curated pass-1 assembly) and copied the order it found
there. So the question moves to pass 1. For each join this prints:
  1. where the two flanking pass-2 contigs align on the RagTag reference (its PAF);
  2. what the curated pass-1 assembly has at that point (final pass-1 AGP):
     one unitig running through it (= hifiasm made the join), or two pieces with a gap;
  3. whether those two pieces were already neighbours in HapHiC's automatic pass-1
     build (= HapHiC made the join) or only after curation (= made by hand).

Usage: join_provenance.py RAGTAG.asm.paf PASS1_FINAL.agp PASS1_AUTO.agp OUT.txt
"""
import collections
import sys

PAF, FINAL, AUTO, OUT = sys.argv[1:5]
CONTIGS = ["h1tg000038l", "h1tg000253l", "h1tg000073l", "h1tg000221l", "h1tg000136l", "h1tg000241l",
           "h2tg000089l", "h2tg000001l"]
JOINS = [  # label, contig before the join, contig after it (order along the final chromosome)
    ("chr1_hap1 L1|P (final 262.930 Mb)", "h1tg000253l", "h1tg000073l"),
    ("chr2_hap2 left gap (final 208.928 Mb)", "h1tg000221l", "h1tg000136l"),
    ("chr2_hap2 right gap (final 225.746 Mb)", "h1tg000136l", "h1tg000241l"),
    ("chr2_hap2 L2|Q switch inside contig h1tg000136l (final ~215.0 Mb = 6.07 Mb into the contig)", "h1tg000136l", 6.07e6),
    ("chr4_hap1 0-85|rest (final 85.526 Mb; pollen: real)", "h2tg000089l", "h2tg000001l"),
]
MIN_BLOCK = 10000
base = lambda n: n.split(":::")[0]


def read_agp(path):
    objs = collections.defaultdict(list)
    for l in open(path):
        if l.startswith("#") or not l.strip():
            continue
        f = l.rstrip("\n").split("\t")
        objs[f[0]].append(f)
    return objs


blocks = collections.defaultdict(list)      # contig -> [(qs, qe, strand, target, ts, te)]
qlen = {}
for l in open(PAF):
    f = l.split("\t")
    if f[0] in CONTIGS and int(f[3]) - int(f[2]) >= MIN_BLOCK:
        qlen[f[0]] = int(f[1])
        blocks[f[0]].append((int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8])))
out = ["WHICH STEP MADE EACH JOIN", "  RagTag PAF: %s" % PAF, "  pass-1 curated AGP: %s" % FINAL,
       "  pass-1 automatic (HapHiC) AGP: %s" % AUTO, "", "1. where the pass-2 contigs sit on the RagTag reference (blocks >= 10 kb)"]
span = {}
for c in CONTIGS:
    bs = sorted(blocks.get(c, []))
    if not bs:
        out.append("  %-12s no alignment blocks in the PAF" % c)
        continue
    by_t = collections.Counter()
    for qs, qe, st, t, ts, te in bs:
        by_t[t] += qe - qs
    t0 = by_t.most_common(1)[0][0]
    on = [b for b in bs if b[3] == t0]
    strand = collections.Counter(b[2] for b in on).most_common(1)[0][0]
    span[c] = (t0, min(b[4] for b in on), max(b[5] for b in on), strand)
    out.append("  %-12s %6.2f Mb  -> %s %.3f-%.3f Mb (%s), %.0f%% of its aligned bp; others: %s" % (
        c, qlen[c] / 1e6, t0, span[c][1] / 1e6, span[c][2] / 1e6, strand, 100 * by_t[t0] / sum(by_t.values()),
        ", ".join("%s %.1f Mb" % (t, n / 1e6) for t, n in by_t.most_common()[1:4]) or "none"))
    if len(by_t) > 1 or c == "h1tg000136l":
        out.append("      along the contig: " + "; ".join(
            "q%.2f-%.2f->%s:%.2f-%.2f%s" % (qs / 1e6, qe / 1e6, t, ts / 1e6, te / 1e6, st) for qs, qe, st, t, ts, te in bs[:12]))

fin, auto = read_agp(FINAL), read_agp(AUTO)
where_auto = {}
for o, fs in auto.items():
    ws = [f for f in fs if f[4] not in ("N", "U")]
    for k, f in enumerate(ws):
        where_auto.setdefault(base(f[5]), (o, k, [base(w[5]) for w in ws]))
out += ["", "2-3. at each join: what the curated pass-1 assembly has there, and was it already in HapHiC's automatic build"]
for label, x, y in JOINS:
    out.append("  " + label)
    if not isinstance(y, str):                 # a point inside one contig: carry it onto the reference
        hit = [bk for bk in blocks.get(x, []) if bk[0] <= y <= bk[1]]
        if not hit:
            out.append("     %s has no alignment block covering %.2f Mb into it" % (x, y / 1e6))
            continue
        qs, qe, st, tx, ts, te = hit[0]
        jp = ts + (y - qs) if st == "+" else te - (y - qs)
        out.append("     on %s at ~%.3f Mb" % (tx, jp / 1e6))
        if tx not in fin:
            out.append("     %s is not an object in the pass-1 curated AGP" % tx)
            continue
        for f in [f for f in fin[tx] if int(f[2]) >= jp - 1e6 and int(f[1]) <= jp + 1e6]:
            tag = "gap %s bp %s" % (f[5], " ".join(f[6:9])) if f[4] in ("N", "U") else "%s %s (part %s)" % (f[5], f[8], f[3])
            out.append("       %9.3f-%9.3f  %s%s" % (int(f[1]) / 1e6, int(f[2]) / 1e6, tag,
                                                    "   <-- the switch" if int(f[1]) <= jp <= int(f[2]) else ""))
        continue
    if x not in span or y not in span:
        out.append("     a flanking contig has no placement in the PAF")
        continue
    (tx, xs, xe, _), (ty, ys, ye, _) = span[x], span[y]
    if tx != ty:
        out.append("     the two contigs sit on different reference scaffolds (%s, %s): RagTag did not join them" % (tx, ty))
        continue
    lo, hi = sorted((xe, ys)) if abs(xe - ys) <= abs(ye - xs) else sorted((ye, xs))
    jp = (lo + hi) / 2
    out.append("     on %s, between %.3f and %.3f Mb (join point ~%.3f Mb)" % (tx, lo / 1e6, hi / 1e6, jp / 1e6))
    if tx not in fin:
        out.append("     %s is not an object in the pass-1 curated AGP (objects: %s)" % (tx, ", ".join(list(fin)[:6])))
        continue
    near = [f for f in fin[tx] if int(f[2]) >= jp - 1e6 and int(f[1]) <= jp + 1e6]
    for f in near:
        tag = "gap %s bp %s" % (f[5], " ".join(f[6:9])) if f[4] in ("N", "U") else "%s %s (part %s)" % (f[5], f[8], f[3])
        mark = "   <-- join point" if int(f[1]) <= jp <= int(f[2]) else ""
        out.append("       %9.3f-%9.3f  %s%s" % (int(f[1]) / 1e6, int(f[2]) / 1e6, tag, mark))
    inside = [f for f in near if f[4] not in ("N", "U") and int(f[1]) + 50000 <= jp <= int(f[2]) - 50000]
    if inside:
        out.append("     -> the join point lies inside one piece, %s: the sequence itself runs across it (hifiasm, or a curation"
                   " cut/join inside that unitig if it carries a :::fragment suffix)" % inside[0][5])
        continue
    ws = [f for f in near if f[4] not in ("N", "U")]
    left = [f for f in ws if int(f[2]) <= jp + 50000]
    right = [f for f in ws if int(f[1]) >= jp - 50000]
    if not left or not right:
        out.append("     -> could not find the pieces on both sides of the join point")
        continue
    a, b = base(left[-1][5]), base(right[0][5])
    out.append("     pieces either side: %s | %s" % (a, b))
    wa, wb = where_auto.get(a), where_auto.get(b)
    if not wa or not wb:
        out.append("     -> %s not found in the automatic build" % (a if not wa else b))
        continue
    adj = wa[0] == wb[0] and abs(wa[1] - wb[1]) == 1
    out.append("     automatic build: %s is in %s (neighbours %s), %s is in %s (neighbours %s)" % (
        a, wa[0], ", ".join(wa[2][max(0, wa[1] - 1):wa[1] + 2]), b, wb[0], ", ".join(wb[2][max(0, wb[1] - 1):wb[1] + 2])))
    out.append("     -> %s" % ("ALREADY NEIGHBOURS in HapHiC's automatic build: HapHiC made this join" if adj else
                              "NOT neighbours in the automatic build: the join was made during curation"))
open(OUT, "w").write("\n".join(out) + "\n")
print("\n".join(out))
