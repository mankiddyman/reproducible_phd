#!/usr/bin/env python3
# Supermatrix per NEPENTHES INTERVAL. Taxa = species:chromosome present in that
# interval. All taxa share the same anchor genes -> no disjoint-data artifact.
import os, csv, collections
W, MINOCC, MINTAX = 25, 0.30, 5
GSD = "/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/results"

ordm = {}
with open(os.path.join(GSD, "combBed.txt")) as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for x in r:
        if x["genome"] == "Nepenthes_gracilis":
            ordm[x["id"]] = int(x["ord"])

rows = list(csv.DictReader(open("super/tip_block.tsv"), delimiter="\t"))
byreg = collections.defaultdict(set)
tipinfo = collections.defaultdict(list)
for r in rows:
    byreg[r["nep_chr"]].add(r["nep_gene"])
    if r["block"] != "NEP":
        ab = {"Dionaea_muscipula":"Dio","Drosera_regia":"reg","Drosera_binata":"bin"}
        tx = ab.get(r["genome"], r["genome"][:3]) + "_" + \
             r["chr"].replace("_hap1","").replace("_collapsed","").replace("chr","c")
        tipinfo[r["nep_gene"]].append((r["tip"], tx))
    else:
        tipinfo[r["nep_gene"]].append((r["tip"], "NEP"))

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

nwin = 0
for reg, genes in sorted(byreg.items()):
    gs = sorted([g for g in genes if g in ordm], key=lambda g: ordm[g])
    for wi in range(0, len(gs), W):
        chunk = gs[wi:wi+W]
        if len(chunk) < W // 2: continue
        alns, taxa = [], set()
        for g in chunk:
            f = f"wgd/aln/{g}.aln"
            if not os.path.exists(f): continue
            fa = read_fa(f)
            if not fa: continue
            best = {}
            for tip, tx in tipinfo[g]:
                if tip not in fa: continue
                gap = fa[tip].count("-")
                if tx not in best or gap < best[tx][0]: best[tx] = (gap, fa[tip])
            if len(best) < MINTAX: continue
            alns.append((len(next(iter(fa.values()))), best))
            taxa |= set(best)
        if len(alns) < 8: continue
        tot = sum(L for L, _ in alns)
        seqs = {t: [] for t in taxa}
        for L, best in alns:
            for t in taxa:
                seqs[t].append(best[t][1] if t in best else "-"*L)
        keep = [t for t in taxa
                if 1 - "".join(seqs[t]).count("-")/tot >= MINOCC or t == "NEP"]
        if len(keep) < MINTAX: continue
        name = f"{reg}_w{wi//W:02d}"
        with open(f"win25/aln/{name}.fa", "w") as fh:
            for t in sorted(keep): fh.write(f">{t}\n{''.join(seqs[t])}\n")
        nwin += 1
        print(f"{name}: {len(keep)} taxa, {len(alns)} genes, {tot} cols, "
              f"[{' '.join(sorted(keep))}]")
print(f"\n{nwin} windows written")
