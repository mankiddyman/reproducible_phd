#!/usr/bin/env python3
# Per ancestral region: concatenate all per-anchor alignments into a supermatrix
# where each TAXON IS A SYNTENIC BLOCK. One tip per block, hundreds of genes deep.
import os, csv, collections

ALN, OUT = "wgd/aln", "super"
rows = list(csv.DictReader(open("super/tip_block.tsv"), delimiter="\t"))

byreg = collections.defaultdict(list)
for r in rows:
    byreg[r["nep_chr"]].append(r)

def read_fa(p):
    d, n, s = {}, None, []
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            if n: d[n] = "".join(s)
            n, s = line[1:].split()[0], []
        else:
            s.append(line)
    if n: d[n] = "".join(s)
    return d

for reg, rr in sorted(byreg.items()):
    anchors = sorted({r["nep_gene"] for r in rr})
    blocks  = sorted({r["block"] for r in rr})
    seqs = {b: [] for b in blocks}
    parts, pos, used = [], 1, 0
    for a in anchors:
        f = os.path.join(ALN, a + ".aln")
        if not os.path.exists(f): continue
        fa = read_fa(f)
        if not fa: continue
        L = len(next(iter(fa.values())))
        m = {r["tip"]: r["block"] for r in rr if r["nep_gene"] == a}
        best = {}
        for tip, sq in fa.items():
            b = m.get(tip)
            if b is None: continue
            g = sq.count("-")
            if b not in best or g < best[b][0]: best[b] = (g, sq)
        if len(best) < 4: continue
        for b in blocks:
            seqs[b].append(best[b][1] if b in best else "-" * L)
        parts.append((a, pos, pos + L - 1)); pos += L; used += 1
    if used < 5: 
        print(f"{reg}: only {used} usable anchors, skipped"); continue
    occ = {b: 1 - "".join(seqs[b]).count("-") / pos for b in blocks}
    keepb = [b for b in blocks if occ[b] >= 0.15 or b == "NEP"]
    with open(f"{OUT}/{reg}.fa", "w") as fh:
        for b in keepb:
            fh.write(f">{b}\n{''.join(seqs[b])}\n")
    with open(f"{OUT}/{reg}.part", "w") as fh:
        fh.write("#nexus\nbegin sets;\n")
        for i, (a, s, e) in enumerate(parts, 1):
            fh.write(f"  charset p{i} = {s}-{e};\n")
        fh.write("end;\n")
    print(f"{reg}: {len(keepb)} blocks, {used} genes, {pos-1} aa columns")
