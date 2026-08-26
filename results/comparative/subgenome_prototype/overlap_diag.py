#!/usr/bin/env python3
# How much data do the supermatrix tips actually share?
import glob, os, itertools
def read_fa(p):
    d, n, s = {}, None, []
    for l in open(p):
        l = l.rstrip()
        if l.startswith(">"):
            if n: d[n] = "".join(s)
            n, s = l[1:].split()[0], []
        else: s.append(l)
    if n: d[n] = "".join(s)
    return d
print(f"{'region':<12}{'tips':>5}{'cols':>8}{'medOverlap':>12}{'pairs<100':>11}")
for f in sorted(glob.glob("super/chr*_dom.fa")):
    fa = read_fa(f); ks = list(fa)
    L = len(fa[ks[0]])
    ov, bad = [], 0
    for a, b in itertools.combinations(ks, 2):
        A, B = fa[a], fa[b]
        n = sum(1 for i in range(L) if A[i] != "-" and B[i] != "-")
        ov.append(n)
        if n < 100: bad += 1
    ov.sort()
    print(f"{os.path.basename(f)[:-3]:<12}{len(ks):>5}{L:>8}"
          f"{ov[len(ov)//2]:>12}{bad}/{len(ov):>10}")
