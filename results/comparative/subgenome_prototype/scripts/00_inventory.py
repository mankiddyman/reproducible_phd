#!/usr/bin/env python3
"""Stage 00 — inventory. READ ONLY. Fails loudly rather than letting a bad join
propagate into the trees."""
import os, sys, glob, csv, shutil, statistics as st
from collections import Counter, defaultdict

BASE = os.environ.get("SUBG_BASE", os.getcwd())
CODON_DIR = os.path.join(BASE, "wgd7", "codon_gap")
META = os.path.join(BASE, "wgd7", "tip_meta.tsv")
GENOMES = ["Dionaea_muscipula","Drosera_binata","Drosera_capensis","Drosera_paradoxa",
           "Drosera_regia","Drosera_scorpioides","Nepenthes_gracilis"]
NEP = "Nepenthes_gracilis"
fail, warn = [], []

def rule(t): print(f"\n{'='*68}\n{t}\n{'='*68}")

def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None: yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else: buf.append(line.strip())
    if name is not None: yield name, "".join(buf)

rule("1. tools")
for tool in ["iqtree2","iqtree","parallel","mafft","seqkit","Rscript"]:
    p = shutil.which(tool); print(f"  {'OK ' if p else 'MISSING'}  {tool:10s} {p or ''}")
if not (shutil.which("iqtree2") or shutil.which("iqtree")): fail.append("no iqtree on PATH")
if not shutil.which("parallel"): fail.append("no GNU parallel on PATH")

rule("2. codon alignments (wgd7/codon_gap)")
fnas = sorted(glob.glob(os.path.join(CODON_DIR, "*.fna")))
print(f"  files: {len(fnas)}")
if not fnas:
    print(f"  !! nothing at {CODON_DIR}"); print("\nABORT"); sys.exit(1)

lengths, ntips, bad_mod3, ragged, nep_count = [], [], [], [], Counter()
aln_tips = {}
for f in fnas:
    anchor = os.path.basename(f)[:-4]
    recs = list(read_fasta(f))
    if not recs: ragged.append((anchor,"empty")); continue
    L = {len(s) for _, s in recs}
    if len(L) > 1: ragged.append((anchor, f"unequal {sorted(L)[:3]}")); continue
    L = L.pop(); lengths.append(L); ntips.append(len(recs))
    if L % 3: bad_mod3.append((anchor, L))
    names = [n for n, _ in recs]; aln_tips[anchor] = names
    nep_count[sum(1 for n in names if n.startswith(NEP))] += 1

print(f"  length (nt): min {min(lengths)}  median {int(st.median(lengths))}  max {max(lengths)}")
print(f"  codons:      min {min(lengths)//3}  median {int(st.median(lengths))//3}  max {max(lengths)//3}")
print(f"  tips/aln:    min {min(ntips)}  median {int(st.median(ntips))}  max {max(ntips)}  total {sum(ntips)}")
print(f"  < 100 codons: {sum(1 for L in lengths if L//3 < 100)}")
print(f"  < 4 tips (untreeable): {sum(1 for n in ntips if n < 4)}")
if bad_mod3:
    print(f"  !! {len(bad_mod3)} not divisible by 3, e.g. {bad_mod3[:3]}")
    fail.append("length not divisible by 3 — codon partitioning invalid")
else: print("  OK  all lengths divisible by 3")
if ragged:
    print(f"  !! {len(ragged)} ragged/empty, e.g. {ragged[:3]}"); fail.append("ragged alignments")
print(f"  Nepenthes tips per alignment: {dict(sorted(nep_count.items()))}")
if nep_count.get(1,0) != len(aln_tips):
    fail.append("not every alignment has exactly one Nepenthes tip — rooting rule breaks")
else: print("  OK  exactly one Nepenthes tip everywhere")

rule("3. tip naming")
for s in aln_tips[sorted(aln_tips)[0]][:3]: print(f"    {s}")
sep, unprefixed = Counter(), []
for a, names in aln_tips.items():
    for n in names:
        sep["@" if "@" in n else "no-@"] += 1
        if not any(n.startswith(g) for g in GENOMES): unprefixed.append((a,n))
print(f"  separator counts: {dict(sep)}")
if unprefixed:
    print(f"  !! {len(unprefixed)} tips lack a known genome prefix, e.g. {unprefixed[:3]}")
    fail.append("unrecognised genome prefix — species parsing will fail")
else: print("  OK  every tip starts with one of the 7 genome names")

rule("4. tip_meta.tsv")
if not os.path.exists(META):
    print(f"  !! missing: {META}"); fail.append("tip_meta.tsv missing")
else:
    with open(META) as fh:
        rdr = csv.DictReader(fh, delimiter="\t"); cols = rdr.fieldnames or []; rows = list(rdr)
    print(f"  rows: {len(rows)}\n  columns: {cols}")
    for c in ["nep_gene","genome"]:
        if c not in cols: fail.append(f"tip_meta missing column {c}")
    chrom = [c for c in cols if c.lower() in
             ("chrom","chr","chromosome","block","tip_chr","chr_label","seqid")]
    print(f"  chromosome-like columns: {chrom or 'NONE FOUND'}")
    if not chrom: warn.append("no chromosome column — co-occurrence step needs one")
    tipcol = next((c for c in cols if c.lower() in ("tip","tip_label","id","gene","seq")), None)
    print(f"  tip-label column guess: {tipcol}")
    if tipcol:
        meta_tips = defaultdict(set)
        for r in rows: meta_tips[r["nep_gene"]].add(r[tipcol])
        common = set(aln_tips) & set(meta_tips)
        print(f"  anchors aln {len(aln_tips)}  meta {len(meta_tips)}  shared {len(common)}")
        miss = set(aln_tips) - set(meta_tips)
        if miss:
            print(f"  !! {len(miss)} anchors absent from meta, e.g. {sorted(miss)[:3]}")
            warn.append("anchors missing from tip_meta")
        checked = matched = 0
        for a in sorted(common)[:100]:
            for t in aln_tips[a]:
                checked += 1; matched += t in meta_tips[a]
        rate = matched/checked if checked else 0
        print(f"  tip-level join, first 100 anchors: {matched}/{checked} = {rate:.3f}")
        if rate < 0.99: fail.append(f"tip join rate {rate:.3f} — names do not match tip_meta")
        else: print("  OK  tips join to meta")
    else: warn.append("could not guess tip-label column — inspect header manually")

rule("5. peptide alignments (AA set)")
found = []
for pat in ["wgd7/aln/*","wgd7/pep/*","wgd7/faa/*","wgd7/aa/*","wgd7/codon/*","wgd/aln/*"]:
    hits = glob.glob(os.path.join(BASE, pat))
    if hits: found.append((pat, len(hits), os.path.basename(hits[0])))
if found:
    for pat,n,ex in found: print(f"  {pat:20s} {n:5d} files   e.g. {ex}")
    print("  -> set AA_DIR to whichever holds ALIGNED peptides")
else:
    print("  none found"); warn.append("peptide alignments not located — AA needs mafft rebuild")

rule("6. anchor -> ancestral region")
ex = sorted(aln_tips)[0]; print(f"  example anchor id: {ex}")
if "_dom" in ex or ex.startswith("chr"):
    regions = Counter(a.split("-")[0] for a in aln_tips)
    print(f"  parsed from anchor name: {len(regions)} regions")
    for r,n in sorted(regions.items()): print(f"    {r:16s} {n}")
else:
    print("  !! not parseable from the anchor id"); warn.append("need an anchor->region map")

rule("VERDICT")
if fail:
    print("BLOCKING:"); [print(f"  X  {f}") for f in fail]
if warn:
    print("WARNINGS (non-blocking):"); [print(f"  !  {w}") for w in warn]
print("\n  No blocking problems. Proceed to 01_prep.sh" if not fail
      else "\n  Fix blocking problems before building.")
sys.exit(1 if fail else 0)
