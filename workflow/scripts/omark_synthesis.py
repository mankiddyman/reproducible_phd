#!/usr/bin/env python3
"""Synthesise OMArk into two static tables. Table 1: master per (species,tool),
Viridiplantae-based, completeness + consistency + U-decomposition (% of Viridi-U,
the three classes sum to 100% of U). Table 2: recurring nonplant clades across runs."""
import glob, re, csv, os, collections
NONPLANT_KINGDOMS = {"Metazoa","Fungi","Bacteria","OtherEukaryote","Archaea"}

def parse_sum(path):
    t=open(path).read()
    npro=re.search(r"there are (\d+) proteins",t)
    cl=re.search(r"^S:([\d.]+)%,D:([\d.]+)%\[.*?\],M:([\d.]+)%",t,re.M)
    co=re.search(r"^A:([\d.]+)%\[.*?\],I:([\d.]+)%\[.*?\],C:([\d.]+)%\[.*?\],U:([\d.]+)%",t,re.M)
    mg=re.search(r"Magnoliopsida\s+\d+\s+\d+\s+([\d.]+)%",t)
    return {"n_prot":int(npro.group(1)) if npro else 0,
        "S":float(cl.group(1)),"D":float(cl.group(2)),"M":float(cl.group(3)),
        "A":float(co.group(1)),"I":float(co.group(2)),"C":float(co.group(3)),"U":float(co.group(4)),
        "magno":float(mg.group(1)) if mg else 0.0}

def parse_rollup(path):
    king,clade={},{}; n_unknown=0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("species\t"): continue
            f=line.rstrip("\n").split("\t")
            if len(f)<6: continue
            nu,level,cat,count=f[2],f[3],f[4],f[5]
            n_unknown=int(nu); count=int(count)
            if level=="kingdom": king[cat]=count
            elif level=="classification" and cat.startswith("nonplant:"):
                clade[cat.split("nonplant:",1)[1]]=count
    return n_unknown,king,clade

rows=[]; clade_runs=collections.Counter(); clade_total=collections.Counter()
for summ in sorted(glob.glob("results/*/annotation/*/omark/viridiplantae/*.viridiplantae.sum")):
    m=re.search(r"results/([^/]+)/annotation/([^/]+)/omark/",summ)
    sp,tool=m.group(1),m.group(2)
    s=parse_sum(summ)
    rollup=f"results/{sp}/annotation/{tool}/omark/{sp}.{tool}.rollup.tsv"
    if not os.path.exists(rollup): continue
    n_unknown,king,clade=parse_rollup(rollup)
    plant_orphan=king.get("Viridiplantae",0); unplaceable=king.get("unplaceable",0)
    nonplant=sum(king.get(k,0) for k in NONPLANT_KINGDOMS)
    denom=n_unknown if n_unknown else 1
    pc=lambda x:f"{100.0*x/denom:.1f}%,{x}"
    rows.append({"species":sp.replace("Drosera_","D_"),"tool":tool,"n_prot":s["n_prot"],
        "S":s["S"],"D":s["D"],"M":s["M"],"A":s["A"],"I":s["I"],"C":s["C"],"U":s["U"],"magno":s["magno"],
        "U_n":n_unknown,"nonplant":pc(nonplant),"plant_orphan":pc(plant_orphan),
        "unplaceable":pc(unplaceable),"nonplant_n":nonplant})
    for c,n in clade.items(): clade_runs[c]+=1; clade_total[c]+=n
rows.sort(key=lambda r:(r["species"],r["tool"]))

cols=[("species",10),("tool",16),("n_prot",7),("S",5),("D",5),("M",5),("A",5),("I",5),
      ("C",4),("U",6),("magno",6),("U_n",6),("nonplant",11),("plant_orphan",12),
      ("unplaceable",13),("nonplant_n",7)]
cell=lambda v:f"{v:.1f}" if isinstance(v,float) else str(v)
print("="*60+"\nTABLE 1 - OMArk master (Viridiplantae)\n"+"="*60)
print("  ".join(h.ljust(w) for h,w in cols)); print("  ".join("-"*w for _,w in cols))
for r in rows: print("  ".join(cell(r[h]).ljust(w) for h,w in cols))
with open("results/omark_master.csv","w",newline="") as fh:
    cw=csv.writer(fh); cw.writerow([h for h,_ in cols])
    for r in rows: cw.writerow([r[h] for h,_ in cols])

print("\n"+"="*60+"\nTABLE 2 - recurring nonplant clades (across all 8 runs)\n"+"="*60)
print(f"{'clade':35} {'run_count':>9} {'total_genes':>12}"); print("-"*58)
for c in sorted(clade_total,key=lambda c:(-clade_total[c],-clade_runs[c]))[:30]:
    print(f"{c[:35]:35} {clade_runs[c]:>9} {clade_total[c]:>12}")
with open("results/omark_nonplant_clades.csv","w",newline="") as fh:
    cw=csv.writer(fh); cw.writerow(["clade","run_count","total_genes"])
    for c in sorted(clade_total,key=lambda c:(-clade_total[c],-clade_runs[c])):
        cw.writerow([c,clade_runs[c],clade_total[c]])
print("\n-> results/omark_master.csv  results/omark_nonplant_clades.csv")
