#!/usr/bin/env python3
"""hifi_junction_check.py -- do HiFi reads run straight through each suspect join?

A HiFi read is one physical molecule of 10-25 kb. Where a join is real, reads align
straight across it. Where a contig is chimeric, reads stop at the join (clipped) and
the clipped part aligns where that molecule really continues. At a scaffold gap no
read can cross, but the clipped parts at the two gap edges still say where each end
really continues.

Reads: the CO_smk HiFi alignment to the crossover reference (chr1_hap1, chr2_hap2,
chr3-6_hap1; minimap2 map-hifi, both haplotypes' reads on one copy).

Usage: hifi_junction_check.py HIFI.bam OUT.txt
"""
import collections
import re
import sys

import pysam

BAM, OUT = sys.argv[1], sys.argv[2]
WINDOWS = [  # chrom, start Mb, end Mb, label
    ("chr1_hap1", 255.0, 272.0, "chr1_hap1, L1|P join (RagTag gap at 262.930)"),
    ("chr2_hap2", 203.0, 231.0, "chr2_hap2, L2|Q join (switch near 215, inside contig h1tg000136l 208.928-225.746)"),
    ("chr4_hap1", 80.0, 92.0, "chr4_hap1, 0-85|rest join (RagTag gap at 85.526; pollen: real)"),
    ("chr3_hap1", 140.0, 152.0, "control: an ordinary stretch of chr3_hap1"),
]
SPAN_STEP = 0.5          # Mb between spanning-read counts
FLANK = 5000             # a spanning read covers the point +- 5 kb
MINCLIP, MINMAPQ = 2000, 20
CLUSTER_KB = 20


def aligned_len(cigar):
    return sum(int(n) for n, op in re.findall(r"(\d+)([MIX=])", cigar))


out = ["HOW TO READ THIS",
       "  At chr1_hap1 ~262.9 Mb (end of L1 / start of P, a scaffold gap) reads must stop; what matters is where",
       "  their clipped part lands. L1-end reads landing on chr2_hap2 near 215-226 Mb: L1 continues into Q",
       "  (the pollen's reading). Landing on chr1_hap1 near 263-265 Mb: L1 continues into P (the scaffold as built).",
       "  About half each: a real heterozygous translocation (one homolog each way).",
       "  Inside contig h1tg000136l (chr2_hap2 208.9-225.7 Mb): if the spanning count drops to ~0 at one point and",
       "  reads stop there, the contig itself joins two pieces that are not joined in the plant (a hifiasm chimera).", ""]
bam = pysam.AlignmentFile(BAM)
for c, a, b, label in WINDOWS:
    out.append("=" * 100)
    out.append("%s   (%s:%g-%g Mb)" % (label, c, a, b))
    # 1) spanning reads: one molecule covering the point and 5 kb either side
    row = []
    x = a
    while x <= b + 1e-9:
        p = int(x * 1e6)
        n = sum(1 for r in bam.fetch(c, p - FLANK, p + FLANK)
                if not (r.is_secondary or r.is_supplementary or r.is_unmapped) and r.mapping_quality >= MINMAPQ
                and r.reference_start <= p - FLANK and r.reference_end >= p + FLANK)
        row.append((x, n))
        x += SPAN_STEP
    out.append("  reads spanning each point (+-5 kb, MAPQ >= %d), every %g Mb:" % (MINMAPQ, SPAN_STEP))
    for k in range(0, len(row), 8):
        out.append("    " + "  ".join("%7.1f:%-4d" % (x, n) for x, n in row[k:k + 8]))
    # 2) reads clipped by >= 2 kb, and where their clipped part aligns
    clusters = collections.defaultdict(lambda: collections.Counter())
    for r in bam.fetch(c, int(a * 1e6), int(b * 1e6)):
        if r.is_secondary or r.is_unmapped or r.mapping_quality < MINMAPQ or not r.cigartuples:
            continue
        for side, (op, ln) in (("left", r.cigartuples[0]), ("right", r.cigartuples[-1])):
            if op not in (4, 5) or ln < MINCLIP:
                continue
            pos = r.reference_start if side == "left" else r.reference_end
            if not (a * 1e6 <= pos <= b * 1e6):
                continue
            dest = "no other alignment"
            if r.has_tag("SA"):
                sa = [s.split(",") for s in r.get_tag("SA").strip(";").split(";") if s]
                sa = [s for s in sa if len(s) >= 6]
                if sa:
                    best = max(sa, key=lambda s: aligned_len(s[3]))
                    dest = "%s:%.1f Mb (MAPQ %s)" % (best[0], int(best[1]) / 1e6,
                                                     "<20" if int(best[4]) < 20 else ">=20")
            clusters[(int(pos / (CLUSTER_KB * 1e3)), side)][dest] += 1
    hot = sorted(((sum(v.values()), k, v) for k, v in clusters.items() if sum(v.values()) >= 5), key=lambda t: -t[0])[:15]
    out.append("  the 15 places where most reads stop (>= 5 reads clipped by >= %d kb), and where their clipped part aligns:" % (
        MINCLIP // 1000))
    shown = 0
    for n, (binx, side), dests in sorted(hot, key=lambda t: t[1][0]):
        shown += 1
        out.append("    %8.3f Mb  %-5s  %3d reads   ->  %s" % (
            binx * CLUSTER_KB / 1e3, side, n, "; ".join("%s x%d" % (d, m) for d, m in dests.most_common(3))))
    if not shown:
        out.append("    none")
open(OUT, "w").write("\n".join(out) + "\n")
print("\n".join(out))
