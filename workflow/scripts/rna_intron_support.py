#!/usr/bin/env python3
"""RNA intron-support: fraction of a tool's predicted introns confirmed by STAR
splice junctions (>=N uniq-mapped reads). precision = of predicted introns, %
STAR-confirmed; recall = of STAR junctions, % predicted. Introns from `exon`
grouped by (chr,Parent), deduped. Exact (chr,start,end) match, strand-ignored."""
import sys, re, argparse
def load_star(sj, min_uniq):
    truth=set()
    for line in open(sj):
        f=line.split("\t")
        if len(f)<7: continue
        if int(f[6])>=min_uniq: truth.add((f[0], int(f[1]), int(f[2])))
    return truth
def gff_introns(gff):
    tx={}
    for line in open(gff):
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        if len(f)<9 or f[2]!="exon": continue
        m=re.search(r"Parent=([^;]+)", f[8])
        if not m: continue
        tx.setdefault((f[0], m.group(1)), []).append((int(f[3]), int(f[4])))
    introns=set()
    for (chrom,_), exons in tx.items():
        if len(exons)<2: continue
        exons.sort()
        for i in range(len(exons)-1):
            istart=exons[i][1]+1; iend=exons[i+1][0]-1
            if iend>=istart: introns.add((chrom, istart, iend))
    return introns
ap=argparse.ArgumentParser()
ap.add_argument("--gff", required=True); ap.add_argument("--sj", required=True)
ap.add_argument("--min-uniq", type=int, default=3)
ap.add_argument("--species", default="?"); ap.add_argument("--tool", default="?")
ap.add_argument("--validate", action="store_true")
a=ap.parse_args()
truth=load_star(a.sj, a.min_uniq); pred=gff_introns(a.gff); conf=pred & truth
prec=100.0*len(conf)/len(pred) if pred else 0.0
rec=100.0*len(conf)/len(truth) if truth else 0.0
if a.validate:
    sys.stderr.write(f"# {a.tool}: n_pred={len(pred)} n_truth={len(truth)} confirmed={len(conf)}\n")
    for x in sorted(pred)[:5]: sys.stderr.write(f"#   {x}  in_truth={x in truth}\n")
print(f"{a.species}\t{a.tool}\t{len(pred)}\t{len(conf)}\t{prec:.1f}\t{rec:.1f}\t{a.min_uniq}\t{len(truth)}")
