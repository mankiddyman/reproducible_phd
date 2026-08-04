#!/usr/bin/env python3
"""Per-gene keep/drop decision from OMArk + the LUCA merge.

KEEP transcript if: Consistent (any P/F) OR Unknown-attributed plant_orphan/unplaceable.
DROP transcript if: Inconsistent OR Contamination OR Unknown-attributed nonplant.
Gene-level rollup: DROP a gene iff ALL its transcripts are DROP (keep if >=1 keeps).

Inputs:
  --ump          viridiplantae .ump (>Consistent_* / >Inconsistent_* / >Contamination_* / >Unknown)
  --attribution  merge attribution.tsv (protein_id ... classification)  [splits the Unknown set]
  --gff          tool GFF3 (to enumerate genes<->transcripts; authoritative gene/tx structure)
Output TSV: gene_id  n_tx  n_keep_tx  decision  reason
"""
import sys, re, argparse, collections

def parse_ump(path):
    """transcript_id -> top-level class in {Consistent,Inconsistent,Contamination,Unknown}"""
    cls={}; cur=None
    for line in open(path):
        line=line.rstrip("\n")
        if line.startswith(">"):
            cur=line[1:].split("_")[0]   # Consistent_Full -> Consistent
        elif line.strip():
            cls[line.strip()]=cur
    return cls

def parse_attr(path):
    """transcript_id -> merge classification (plant_orphan|unplaceable|nonplant:*)"""
    out={}
    with open(path) as fh:
        for i,line in enumerate(fh):
            f=line.rstrip("\n").split("\t")
            if i==0 and f[0]=="protein_id": continue
            if len(f)<7: continue
            out[f[0]]=f[6]
    return out

def gene_tx_map(gff):
    """gene_id -> [transcript_ids] from mRNA/transcript Parent; tx_id -> gene_id too."""
    g2t=collections.defaultdict(list); t2g={}
    for line in open(gff):
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        if len(f)<9 or f[2] not in ("mRNA","transcript"): continue
        mi=re.search(r"ID=([^;]+)",f[8]); mp=re.search(r"Parent=([^;]+)",f[8])
        if not mi or not mp: continue
        tx=mi.group(1); gene=mp.group(1)
        g2t[gene].append(tx); t2g[tx]=gene
    return g2t, t2g

ap=argparse.ArgumentParser()
ap.add_argument("--ump", required=True); ap.add_argument("--attribution", required=True)
ap.add_argument("--gff", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--species", default="?"); ap.add_argument("--tool", default="?")
ap.add_argument("--rna-evidenced", action="store_true",
                help="tool used RNA evidence (e.g. braker_final): keep unscored loci as RNA-supported; else drop")
a=ap.parse_args()

RNA_EVIDENCED=a.rna_evidenced
ump=parse_ump(a.ump)
attr=parse_attr(a.attribution)
g2t,t2g=gene_tx_map(a.gff)

def tx_decision(tx):
    c=ump.get(tx)
    if c in ("Inconsistent","Contamination"): return "drop", c
    if c=="Consistent": return "keep","Consistent"
    if c=="Unknown":
        cl=attr.get(tx,"unplaceable")          # default unplaceable if merge silent
        if cl.startswith("nonplant"): return "drop","nonplant"
        if cl=="plant_orphan": return "keep","plant_orphan"
        return "keep","unplaceable"
    if RNA_EVIDENCED:
        return "keep","unscored_rna"            # RNA-built locus, ORF too short for omamer -> keep (transcribed)
    return "drop","unscored"                    # ab-initio, unscored & too short -> spurious, drop

# reconciliation sanity
n_ump=len(ump); n_attr=len(attr)
n_unknown=sum(1 for v in ump.values() if v=="Unknown")
sys.stderr.write(f"# {a.tool}: ump_tx={n_ump} unknown_in_ump={n_unknown} attribution_rows={n_attr}\n")
if n_unknown!=n_attr:
    sys.stderr.write(f"# WARN: Unknown count ({n_unknown}) != attribution rows ({n_attr})\n")

rows=[]; kc=collections.Counter()
for gene,txs in g2t.items():
    decs=[tx_decision(t) for t in txs]
    n_keep=sum(1 for d,_ in decs if d=="keep")
    decision = "keep" if n_keep>=1 else "drop"
    # reason: dominant keep-reason if kept, else dominant drop-reason
    reasons=[r for d,r in decs if d==decision]
    reason=collections.Counter(reasons).most_common(1)[0][0] if reasons else "NA"
    rows.append((gene,len(txs),n_keep,decision,reason))
    kc[(decision,reason)]+=1

with open(a.out,"w") as fh:
    fh.write("gene_id\tn_tx\tn_keep_tx\tdecision\treason\n")
    for r in rows: fh.write("\t".join(map(str,r))+"\n")

# summary to stderr
sys.stderr.write(f"# genes: {len(rows)}  kept: {sum(1 for r in rows if r[3]=='keep')}  dropped: {sum(1 for r in rows if r[3]=='drop')}\n")
for (d,r),n in sorted(kc.items()): sys.stderr.write(f"#   {d:5} {r:15} {n}\n")
