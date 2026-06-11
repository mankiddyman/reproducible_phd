#!/usr/bin/env python3
"""Attribute Viridiplantae-unknown proteins via the LUCA run.

For each protein OMArk(Viridiplantae) calls 'Unknown', look up where
OMArk(LUCA) placed it (omamer hoglevel) and classify:
  - placed in a clade descending from Viridiplantae (33090) -> plant_orphan  (KEEP: real plant gene the Viridiplantae HOG set under-credited)
  - placed in a non-plant clade                              -> nonplant:<clade> (DROP: non-target. NOT asserted as contamination -- could be contaminant, shared-reference contaminant, or distant homology; indistinguishable by score)
  - no LUCA placement (hoglevel N/A)                         -> unplaceable  (KEEP: novel / Drosera-specific candidate or model noise; no contamination evidence)

Outputs:
  --out-detail   per-protein TSV (id, hoglevel, taxid, kingdom, family_p, overlap, classification)
  --out-rollup   classification + kingdom summary
  --out-droplist .txt of protein IDs to REMOVE (the nonplant set) -> consumed by the proteome/GFF3 refinement rule
  --out-orphans  the plant_orphan set with HOGs (post-seminar curiosity: LUCA-rescued plant genes Viridiplantae missed)
"""
import argparse, sys, collections
from ete3 import NCBITaxa

VIRIDIPLANTAE = 33090
# kingdom anchors, checked in lineage order; 2759 (Eukaryota) is the catch-all last
KINGDOM_ANCHORS = [(33090,"Viridiplantae"),(33208,"Metazoa"),(4751,"Fungi"),
                   (2,"Bacteria"),(2157,"Archaea"),(10239,"Viruses"),(2759,"OtherEukaryote")]

def parse_ump_unknown(path):
    unknown, cur = set(), None
    for line in open(path):
        line = line.rstrip("\n")
        if not line: continue
        if line.startswith(">"): cur = line[1:]; continue
        if cur == "Unknown": unknown.add(line.strip())
    return unknown

def parse_omamer(path):
    """qseqid -> (hoglevel, family_p, overlap). N/A-placed proteins kept as ('N/A',nan,nan)."""
    d = {}
    for line in open(path):
        if line.startswith("!"): continue
        c = line.rstrip("\n").split("\t")
        if not c or c[0] == "qseqid": continue
        hog = c[2] if len(c) > 2 else "N/A"
        try: fam_p = float(c[3])
        except (IndexError, ValueError): fam_p = float("nan")
        try: ov = float(c[10])
        except (IndexError, ValueError): ov = float("nan")
        d[c[0]] = (hog, fam_p, ov)
    return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--viridi-ump", required=True)
    ap.add_argument("--luca-omamer", required=True)
    ap.add_argument("--out-detail", required=True)
    ap.add_argument("--out-rollup", required=True)
    ap.add_argument("--out-droplist", required=True)
    ap.add_argument("--out-orphans", required=True)
    ap.add_argument("--species", default="NA")
    ap.add_argument("--tool", default="NA")
    args = ap.parse_args()

    ncbi = NCBITaxa()
    unknown = parse_ump_unknown(args.viridi_ump)
    luca = parse_omamer(args.luca_omamer)

    clade_cache, king_cache = {}, {}
    def resolve(clade):
        if clade in clade_cache: return clade_cache[clade]
        if clade in ("N/A","",None): r=(None,False); clade_cache[clade]=r; return r
        tr = ncbi.get_name_translator([clade])
        if clade not in tr: r=(None,None)  # unresolved name -> treat as unplaceable
        else:
            tid = tr[clade][0]; lin = set(ncbi.get_lineage(tid))
            r=(tid, VIRIDIPLANTAE in lin)
        clade_cache[clade]=r; return r
    def kingdom(clade):
        if clade in king_cache: return king_cache[clade]
        tr = ncbi.get_name_translator([clade])
        if clade not in tr: k="unplaceable"
        else:
            lin = set(ncbi.get_lineage(tr[clade][0])); k="other"
            for tid,lab in KINGDOM_ANCHORS:
                if tid in lin: k=lab; break
        king_cache[clade]=k; return k

    rows, droplist, orphans = [], [], []
    cls_roll, king_roll = collections.Counter(), collections.Counter()
    for prot in sorted(unknown):
        hog, fam_p, ov = luca.get(prot, ("N/A", float("nan"), float("nan")))
        tid, is_plant = resolve(hog)
        if hog in ("N/A","") or is_plant is None:
            cls, king = "unplaceable", "unplaceable"
        elif is_plant:
            cls, king = "plant_orphan", "Viridiplantae"
            orphans.append((prot, hog, fam_p, ov))
        else:
            cls, king = f"nonplant:{hog}", kingdom(hog)
            droplist.append(prot)
        rows.append((prot, hog, tid if tid else "", king,
                     f"{fam_p:.2f}" if fam_p==fam_p else "NA",
                     f"{ov:.3f}" if ov==ov else "NA", cls))
        cls_roll[cls]+=1; king_roll[king]+=1

    with open(args.out_detail,"w") as fh:
        fh.write("protein_id\tluca_hoglevel\tluca_taxid\tkingdom\tfamily_p\toverlap\tclassification\n")
        for r in rows: fh.write("\t".join(map(str,r))+"\n")

    with open(args.out_rollup,"w") as fh:
        fh.write("# kingdom rollup\n")
        fh.write("species\ttool\tn_unknown\tlevel\tcategory\tcount\n")
        for k,v in king_roll.most_common():
            fh.write(f"{args.species}\t{args.tool}\t{len(unknown)}\tkingdom\t{k}\t{v}\n")
        for c,v in cls_roll.most_common():
            fh.write(f"{args.species}\t{args.tool}\t{len(unknown)}\tclassification\t{c}\t{v}\n")

    with open(args.out_droplist,"w") as fh:
        for p in droplist: fh.write(p+"\n")

    with open(args.out_orphans,"w") as fh:
        fh.write("protein_id\tluca_hoglevel\tfamily_p\toverlap\n")
        for p,h,fp,ov in orphans:
            fh.write(f"{p}\t{h}\t{fp:.2f}\t{ov:.3f}\n")

    print(f"[{args.species} {args.tool}] {len(unknown)} unknown -> "
          f"DROP {len(droplist)} nonplant | KEEP {len(orphans)} plant_orphan + "
          f"{cls_roll['unplaceable']} unplaceable", file=sys.stderr)
    print("  kingdom:", dict(king_roll.most_common()), file=sys.stderr)

if __name__ == "__main__":
    main()
