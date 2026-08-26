#!/usr/bin/env python3
"""Pairwise protein p-distance and Nei-Gojobori dN/dS for each triplet."""
import sys, os, math, glob, csv

CODON = {}
BASES = "TCAG"
AAS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
for i,b1 in enumerate(BASES):
    for j,b2 in enumerate(BASES):
        for k,b3 in enumerate(BASES):
            CODON[b1+b2+b3] = AAS[i*16+j*4+k]

def syn_sites(c):
    """expected synonymous sites in a codon (NG86)"""
    if c not in CODON: return 0.0
    s = 0.0
    for pos in range(3):
        same = 0
        for b in BASES:
            if b == c[pos]: continue
            m = c[:pos]+b+c[pos+1:]
            if CODON.get(m) == CODON[c]: same += 1
        s += same/3.0
    return s

def read_fasta(p):
    d, name = {}, None
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"): name = line[1:].split()[0]; d[name] = []
        elif name: d[name].append(line)
    return {k: "".join(v) for k, v in d.items()}

def ng86(a, b):
    """returns (dN, dS, n_sites_syn, n_sites_nonsyn) or None"""
    S = N = Sd = Nd = 0.0
    for i in range(0, min(len(a), len(b)) - 2, 3):
        ca, cb = a[i:i+3].upper(), b[i:i+3].upper()
        if len(ca) < 3 or len(cb) < 3: break
        if any(x not in BASES for x in ca+cb): continue
        if CODON[ca] == "*" or CODON[cb] == "*": continue
        s = (syn_sites(ca) + syn_sites(cb)) / 2.0
        S += s; N += 3 - s
        if ca == cb: continue
        diff = [p for p in range(3) if ca[p] != cb[p]]
        if len(diff) == 1:
            if CODON[ca] == CODON[cb]: Sd += 1
            else: Nd += 1
        else:                       # multi-hit: average over paths, cheap approx
            syn = sum(1 for p in diff if CODON[ca[:p]+cb[p]+ca[p+1:]] == CODON[ca])
            Sd += syn * len(diff)/max(len(diff),1) / len(diff)
            Nd += (len(diff) - syn) * 1.0
    if S < 10 or N < 10: return None
    ps, pn = Sd/S, Nd/N
    def jc(p):
        if p >= 0.749: return float('nan')
        return -0.75*math.log(1-4*p/3)
    return jc(pn), jc(ps), S, N

def pdist(a, b):
    n = d = 0
    for x, y in zip(a, b):
        if x in "-X*" or y in "-X*": continue
        n += 1
        if x != y: d += 1
    return (d/n, n) if n >= 30 else (float('nan'), n)

OUT = sys.argv[1] if len(sys.argv) > 1 else "work"
idx = {r[0]: r[1:] for r in csv.reader(open(f"{OUT}/aln_index.tsv"), delimiter="\t")}
rows = []
for f in sorted(glob.glob(f"{OUT}/aln/*.pal2nal.fna")):
    hid = os.path.basename(f).replace(".pal2nal.fna", "")
    if hid not in idx: continue
    nep, s5, s9, chr1, chr2, arr = idx[hid]
    nt = read_fasta(f)
    aa = read_fasta(f.replace(".pal2nal.fna", ".aln.faa"))
    if not all(k in nt for k in (nep, s5, s9)): continue
    r = {"hog": hid, "nep": nep, "dio1": s5, "dio2": s9, "chr1": chr1,
         "chr2": chr2, "arrangement": arr,
         "chrpair": chr1.split("_")[0] if chr1.split("_")[0]==chr2.split("_")[0] else "mixed",
         "alnlen_nt": len(nt[nep])}
    for lab, x, y in (("dio1_dio2", s5, s9), ("nep_dio1", nep, s5), ("nep_dio2", nep, s9)):
        p, n = pdist(aa.get(x, ""), aa.get(y, ""))
        r[f"pdist_{lab}"] = p; r[f"naa_{lab}"] = n
        g = ng86(nt[x], nt[y])
        r[f"dN_{lab}"], r[f"dS_{lab}"] = (g[0], g[1]) if g else (float('nan'),)*2
    rows.append(r)

cols = ["hog","nep","dio1","dio2","chr1","chr2","chrpair","arrangement","alnlen_nt"] + \
       [f"{m}_{p}" for p in ("dio1_dio2","nep_dio1","nep_dio2") for m in ("pdist","dN","dS","naa")]
with open("distances.csv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
    for r in rows: w.writerow({c: r.get(c,"") for c in cols})
print(f"wrote distances.csv: {len(rows)} triplets")
