#!/usr/bin/env python3
"""Filter a tool's GFF3 by the per-gene OMArk keep/drop table.
Keeps all features of kept genes (gene + mRNA/exon/CDS/UTR/... children),
drops every feature belonging to dropped genes. Tags each kept gene line with
omark_class=<reason>. Preserves IDs and ordering; does not renumber.

Identifies a feature's owning gene via: gene -> ID; everything else -> walk
Parent chain to the gene. Robust to the AGAT structure (gene>mRNA>exon/CDS/...)."""
import sys, re, argparse, collections

def attr(s, key):
    m=re.search(rf"{key}=([^;]+)", s); return m.group(1) if m else None

ap=argparse.ArgumentParser()
ap.add_argument("--gff", required=True); ap.add_argument("--classes", required=True)
ap.add_argument("--out", required=True)
a=ap.parse_args()

# load decisions: gene_id -> (decision, reason)
dec={}
with open(a.classes) as fh:
    for i,line in enumerate(fh):
        f=line.rstrip("\n").split("\t")
        if i==0 or len(f)<5: continue
        dec[f[0]]=(f[3], f[4])

# pass 1: build ID -> parent and ID -> featuretype, to resolve each feature's gene
parent={}; ftype={}; lines=[]
with open(a.gff) as fh:
    for line in fh:
        if line.startswith("#"):
            lines.append(("comment", line)); continue
        f=line.rstrip("\n").split("\t")
        if len(f)<9:
            lines.append(("comment", line)); continue
        fid=attr(f[8],"ID"); par=attr(f[8],"Parent")
        if fid: parent[fid]=par; ftype[fid]=f[2]
        lines.append(("feat", line, f, fid, par))

def owning_gene(fid, par):
    # walk up Parent chain until a feature whose type is gene (or no parent)
    cur_id, cur_par = fid, par
    seen=0
    while cur_par and seen<10:
        if ftype.get(cur_par)=="gene": return cur_par
        cur_id, cur_par = cur_par, parent.get(cur_par)
        seen+=1
    # if no parent, the feature itself may be the gene
    if ftype.get(fid)=="gene": return fid
    return cur_id  # best effort (top of chain)

kept_g=set(g for g,(d,_) in dec.items() if d=="keep")
n_in=n_out=0
with open(a.out,"w") as out:
    for rec in lines:
        if rec[0]=="comment":
            out.write(rec[1]); continue
        _, raw, f, fid, par = rec
        gene = fid if f[2]=="gene" else owning_gene(fid, par)
        n_in+=1
        if gene not in kept_g:
            continue
        if f[2]=="gene":
            reason=dec[gene][1]
            sep="" if f[8].endswith(";") or f[8]=="" else ";"
            f[8]=f[8].rstrip("\n")+sep+f"omark_class={reason}"
            out.write("\t".join(f)+"\n")
        else:
            out.write(raw)
        n_out+=1

sys.stderr.write(f"# features in: {n_in}  out: {n_out}  genes kept: {len(kept_g)}\n")
