#!/usr/bin/env python3
"""Pairwise Ks with PAML yn00 for a chunk of anchors.
yn00 writes scratch files (2YN.dN, rst, rub) into the CWD, so each chunk runs in its
own directory. Sequence names are shortened because PAML dislikes '@' and '.'."""
import os, sys, subprocess, csv, shutil

ROOT = "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype"
ALN  = os.path.join(ROOT, "wgd7/codon_gap")
chunk_file, out_csv, tmpdir = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(tmpdir, exist_ok=True)

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

def parse_yn(path):
    """Read the YN00 pairwise block. Returns {(i,j): (dN, dS)} with 1-based indices."""
    out, inblock = {}, False
    for line in open(path):
        if line.startswith("seq. seq."): inblock = True; continue
        if inblock:
            f = line.split()
            if len(f) < 11:
                if out and not line.strip(): break
                continue
            try:
                i, j = int(f[0]), int(f[1])
                dN, dS = float(f[7]), float(f[10])
                out[(i, j)] = (dN, dS)
            except (ValueError, IndexError):
                continue
    return out

anchors = [l.strip() for l in open(chunk_file) if l.strip()]
cwd = os.getcwd(); os.chdir(tmpdir)
rows = []
for a in anchors:
    fa = os.path.join(ALN, a + ".fna")
    if not os.path.exists(fa): continue
    d = read_fa(fa)
    names = [k for k in d if len(d[k]) > 0]
    if len(names) < 2: continue
    L = len(d[names[0]])
    if L < 150 or L % 3: continue
    short = {f"s{i+1:04d}": n for i, n in enumerate(names)}
    with open("in.phy", "w") as fh:
        fh.write(f"  {len(names)}  {L}\n")
        for sid, n in short.items():
            fh.write(f"{sid}\n{d[n]}\n")
    with open("yn00.ctl", "w") as fh:
        fh.write("seqfile = in.phy\noutfile = out.txt\nverbose = 0\n"
                 "icode = 0\nweighting = 0\ncommonf3x4 = 0\n")
    try:
        subprocess.run(["/usr/bin/yn00"], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, timeout=300, check=False)
        if not os.path.exists("out.txt"): continue
        res = parse_yn("out.txt")
    except Exception:
        continue
    idx = {i+1: names[i] for i in range(len(names))}
    for (i, j), (dN, dS) in res.items():
        a1, a2 = idx.get(i), idx.get(j)
        if a1 is None or a2 is None: continue
        rows.append([a, a1, a2, a1.split("@")[0], a2.split("@")[0], L//3, dN, dS])
    for f in ("out.txt", "2YN.dN", "2YN.dS", "2YN.t", "rst", "rst1", "rub"):
        if os.path.exists(f): os.remove(f)
os.chdir(cwd)
with open(out_csv, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["anchor","seq1","seq2","sp1","sp2","codons","dN","dS"])
    w.writerows(rows)
print(f"{os.path.basename(chunk_file)}: {len(anchors)} anchors -> {len(rows)} pairs")
