#!/usr/bin/env python3
"""compare_orthofinder_runs.py -- did adding species change how OrthoFinder groups the genes of the species two runs
share? (2026-11: the 11-species Droseraceae + Nepenthes run against the 19-species run that adds holocentric /
monocentric pairs from four other clades.)

Two OrthoFinder 3 results directories are compared on the genes of their shared species (genes in both runs;
unassigned genes count as groups of one):

  A  like for like: old Orthogroups.tsv against the new run's HOGs at the most recent common ancestor (MRCA) of the
     shared species. Both mean "all genes descended from one gene in that ancestor", so the differences are what the
     extra species changed (gene trees, rooting, clustering).
  B  old Orthogroups.tsv against new Orthogroups.tsv. In OrthoFinder 3, Orthogroups.tsv holds the HOGs at the root
     (N0). The new root is older, so genes duplicated between the two roots fall into one group by construction:
     this shows what the 19-species tables do to gene families, not a change in clustering.
  C  orthologue pairs between the shared species (inferred from the gene trees): kept, lost, new.
  plus a lookup of named genes (--gene NAME=SPECIES:GENE) in both runs, with per-species gene counts.

Group-level numbers:
  ARI (adjusted Rand index): agreement over all pairs of genes (same group in both runs, or different groups in
      both), corrected for the agreement expected by chance; 1 = identical, about 0 = no better than random.
  pairs kept: share of the old run's same-group gene pairs still together in the new run (falls when groups split).
  pairs new: share of the new run's same-group pairs that were in different groups before (rises when groups merge).
  fate of each old group of 2+ genes: identical | merged (all its genes in one new group that also holds other
      genes) | split (its genes in 2+ new groups that hold nothing else) | reshuffled (anything else).

usage: python3 compare_orthofinder_runs.py --old <Results_dir> --new <Results_dir> --out <dir>
           [--old-log FILE] [--new-log FILE] [--old-input DIR] [--new-input DIR] [--gene NAME=SPECIES:GENE ...]
Standard library only. Writes summary.txt (also printed), old_group_fates.tsv, orthologue_pairs.tsv, gene_lookups.tsv.
"""
import argparse
import collections
import csv
import hashlib
import itertools
import os
import re
import sys

csv.field_size_limit(sys.maxsize)


# ---------------------------------------------------------------------------------------------------- reading
def read_groups(fn, first_sp_col):
    """OrthoFinder group table -> (species columns, {group id: {species: [genes]}}). Cells hold genes joined by ', '."""
    groups = {}
    with open(fn, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        species = header[first_sp_col:]
        for row in rd:
            if not row:
                continue
            g = {}
            for sp, cell in zip(species, row[first_sp_col:]):
                genes = [x.strip() for x in cell.split(",") if x.strip()]
                if genes:
                    g[sp] = genes
            groups[row[0]] = g
    return species, groups


def run_groups(res_dir):
    """Orthogroups.tsv + Orthogroups_UnassignedGenes.tsv of one run -> (species, gene -> group, group -> size)."""
    species, groups = read_groups(os.path.join(res_dir, "Orthogroups", "Orthogroups.tsv"), 1)
    un = os.path.join(res_dir, "Orthogroups", "Orthogroups_UnassignedGenes.tsv")
    if os.path.exists(un):
        _, ug = read_groups(un, 1)
        groups.update(ug)
    gene2grp, sp_of = {}, {}
    for gid, g in groups.items():
        for sp, genes in g.items():
            for x in genes:
                gene2grp[x] = gid
                sp_of[x] = sp
    return species, groups, gene2grp, sp_of


def parse_newick_clades(nwk):
    """Labelled newick -> ({internal label: frozenset(leaves)}, all leaves). Branch lengths/supports are ignored."""
    nwk = nwk.strip().rstrip(";")
    clades, stack, leaves = {}, [[]], []
    i, n = 0, len(nwk)

    def read_name(i):
        j = i
        while j < n and nwk[j] not in "(),:;":
            j += 1
        name = nwk[i:j].strip()
        if j < n and nwk[j] == ":":           # skip a branch length
            j += 1
            while j < n and nwk[j] not in "(),;":
                j += 1
        return name, j

    while i < n:
        c = nwk[i]
        if c == "(":
            stack.append([])
            i += 1
        elif c == ",":
            i += 1
        elif c == ")":
            kids = stack.pop()
            clade = frozenset().union(*kids) if kids else frozenset()
            label, i = read_name(i + 1)
            if label:
                clades[label] = clade
            stack[-1].append(clade)
        else:
            name, i = read_name(i)
            if name:
                leaves.append(name)
                stack[-1].append(frozenset([name]))
    return clades, leaves


def orthologue_pairs(res_dir, sp1, sp2):
    """Set of (gene of sp1, gene of sp2) orthologue pairs from Orthologues/Orthologues_sp1/sp1__v__sp2.tsv, or None."""
    fn = os.path.join(res_dir, "Orthologues", f"Orthologues_{sp1}", f"{sp1}__v__{sp2}.tsv")
    flip = False
    if not os.path.exists(fn):
        fn = os.path.join(res_dir, "Orthologues", f"Orthologues_{sp2}", f"{sp2}__v__{sp1}.tsv")
        flip = True
        if not os.path.exists(fn):
            return None
    pairs = set()
    with open(fn, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd)
        for row in rd:
            if len(row) < 3:
                continue
            a = [x.strip() for x in row[1].split(",") if x.strip()]
            b = [x.strip() for x in row[2].split(",") if x.strip()]
            for x, y in itertools.product(a, b):
                pairs.add((y, x) if flip else (x, y))
    return pairs


def version_line(fn):
    if not fn or not os.path.exists(fn):
        return "(no log given)"
    with open(fn, errors="replace") as fh:
        for line in fh:
            m = re.search(r"OrthoFinder v[0-9][0-9.]*", re.sub(r"\x1b\[[0-9;]*m", "", line))
            if m:
                return m.group(0)
    return "(no version line found)"


def md5(fn):
    h = hashlib.md5()
    with open(fn, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------------------------------------------- comparing
def comb2(x):
    return x * (x - 1) // 2


def compare_partitions(genes, lab_a, lab_b):
    """Pair-counting comparison of two partitions of the same genes. lab_*: gene -> group label."""
    nij = collections.Counter((lab_a[g], lab_b[g]) for g in genes)
    ai = collections.Counter(lab_a[g] for g in genes)
    bj = collections.Counter(lab_b[g] for g in genes)
    s_ij = sum(comb2(v) for v in nij.values())
    s_a = sum(comb2(v) for v in ai.values())
    s_b = sum(comb2(v) for v in bj.values())
    n2 = comb2(len(genes))
    expected = s_a * s_b / n2 if n2 else 0.0
    denom = 0.5 * (s_a + s_b) - expected
    ari = (s_ij - expected) / denom if denom else 1.0
    kept = s_ij / s_a if s_a else 1.0
    new = 1 - s_ij / s_b if s_b else 0.0
    # fate of each old group of 2+ genes
    members = collections.defaultdict(list)
    for g in genes:
        members[lab_a[g]].append(g)
    fates, fate_rows = collections.Counter(), {}
    fate_genes = collections.Counter()
    for grp, gs in members.items():
        if len(gs) < 2:
            continue
        targets = collections.Counter(lab_b[g] for g in gs)
        if len(targets) == 1:
            t = next(iter(targets))
            fate = "identical" if bj[t] == len(gs) else "merged"
        else:
            fate = "split" if all(bj[t] == c for t, c in targets.items()) else "reshuffled"
        fates[fate] += 1
        fate_genes[fate] += len(gs)
        top, top_n = targets.most_common(1)[0]
        fate_rows[grp] = (len(gs), fate, len(targets), round(top_n / len(gs), 3), top, bj[top])
    # genes alone in the old run (unassigned) that joined a group, and the reverse
    joined = sum(1 for g in genes if ai[lab_a[g]] == 1 and bj[lab_b[g]] > 1)
    left = sum(1 for g in genes if ai[lab_a[g]] > 1 and bj[lab_b[g]] == 1)
    n_a = sum(1 for v in ai.values() if v > 1)
    n_b = sum(1 for v in bj.values() if v > 1)
    return dict(ari=ari, kept=kept, new=new, fates=fates, fate_genes=fate_genes, rows=fate_rows,
                joined=joined, left=left, n_groups_a=n_a, n_groups_b=n_b)


def code(sp):
    """Drosera_aliciae -> Dro.alic"""
    p = sp.split("_")
    return f"{p[0][:3]}.{p[1][:4]}" if len(p) > 1 else sp[:8]


def gene_matcher(q):
    """Matches q as a whole token of a gene id: 'g4881' matches g4881.t1, Drosera_regia_g4881.t1, not g48810.t1."""
    rx = re.compile(r"(^|[_|:])" + re.escape(q) + r"($|[.\-_|])")
    return lambda gene: bool(rx.search(gene))


# ---------------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--old", required=True, help="results dir of the run with fewer species (e.g. .../Results_drosera)")
    ap.add_argument("--new", required=True, help="results dir of the run with more species")
    ap.add_argument("--out", required=True)
    ap.add_argument("--old-log"), ap.add_argument("--new-log")
    ap.add_argument("--old-input"), ap.add_argument("--new-input")
    ap.add_argument("--gene", action="append", default=[], help="NAME=SPECIES:GENE (repeatable)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    lines = []

    def say(s=""):
        lines.append(s)
        print(s, flush=True)

    old_sp, old_groups, old_g2, old_spof = run_groups(a.old)
    new_sp, new_groups, new_g2, new_spof = run_groups(a.new)
    shared = [s for s in old_sp if s in set(new_sp)]
    say("OrthoFinder runs compared")
    say(f"  old: {a.old}")
    say(f"       {len(old_sp)} species, {version_line(a.old_log)}")
    say(f"  new: {a.new}")
    say(f"       {len(new_sp)} species, {version_line(a.new_log)}")
    say(f"  shared species ({len(shared)}): {', '.join(shared)}")
    only_old = [s for s in old_sp if s not in set(new_sp)]
    if only_old:
        say(f"  !! species only in the old run (left out of the comparison): {', '.join(only_old)}")
    say(f"  added in the new run ({len(new_sp) - len(shared)}): {', '.join(s for s in new_sp if s not in set(old_sp))}")

    # ---- same genes and proteomes for the shared species?
    say("")
    say("Same input for the shared species?")
    old_by_sp, new_by_sp = collections.defaultdict(set), collections.defaultdict(set)
    for g, s in old_spof.items():
        old_by_sp[s].add(g)
    for g, s in new_spof.items():
        new_by_sp[s].add(g)
    universe = []
    diff_sp = []
    for s in shared:
        common = old_by_sp[s] & new_by_sp[s]
        universe.extend(sorted(common))
        if len(common) != len(old_by_sp[s]) or len(common) != len(new_by_sp[s]):
            diff_sp.append(f"{s} (old {len(old_by_sp[s])}, new {len(new_by_sp[s])}, both {len(common)})")
    if diff_sp:
        say(f"  !! gene sets differ for {len(diff_sp)} species: {'; '.join(diff_sp)}")
    else:
        say(f"  gene ids: identical for all {len(shared)} shared species ({len(universe)} genes)")
    if a.old_input and a.new_input:
        changed = [s for s in shared
                   if os.path.exists(os.path.join(a.old_input, s + ".fa"))
                   and os.path.exists(os.path.join(a.new_input, s + ".fa"))
                   and md5(os.path.join(a.old_input, s + ".fa")) != md5(os.path.join(a.new_input, s + ".fa"))]
        say(f"  input proteomes: {'identical files' if not changed else '!! differ for ' + ', '.join(changed)}")
    universe_set = set(universe)

    # ---- the node of the new species tree matching the old run's root
    clades, leaves = parse_newick_clades(open(os.path.join(a.new, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")).read())
    want = frozenset(shared)
    cands = [(len(c), lab) for lab, c in clades.items() if want <= c]
    node = min(cands)[1]
    extra = sorted(clades[node] - want)
    say("")
    say("Node of the new species tree that matches the old run's root (MRCA of the shared species)")
    say(f"  {node}: {len(clades[node])} species" + (f"; !! it also holds {', '.join(extra)} (the shared species are not one clade)" if extra else ", exactly the shared species"))
    say(f"  new root N0: {len(clades.get('N0', leaves))} species")

    # labels: old groups; new HOGs at the node (genes not in any HOG there: groups of one); new root groups
    lab_old = {g: old_g2[g] for g in universe}
    lab_root = {g: new_g2[g] for g in universe}
    hog_fn = os.path.join(a.new, "Phylogenetic_Hierarchical_Orthogroups", f"{node}.tsv")
    if node == "N0" or not os.path.exists(hog_fn):
        if node != "N0":
            say(f"  !! {hog_fn} not found: using the root groups for both comparisons")
        lab_node = dict(lab_root)
        node_groups = new_groups
    else:
        _, hogs = read_groups(hog_fn, 3)
        node_groups = hogs
        g2hog = {}
        for hid, g in hogs.items():
            for genes in g.values():
                for x in genes:
                    g2hog[x] = hid
        lab_node = {g: g2hog.get(g, ("alone", g)) for g in universe}
        missing = sum(1 for g in universe if g not in g2hog)
        say(f"  HOGs at {node}: {len(hogs)}; genes of the shared species in none of them (counted as groups of one): {missing}")

    # ---- A and B
    res = {}
    for key, title, lab in (("A", f"A  like for like: old groups vs new HOGs at {node}", lab_node),
                            ("B", "B  old groups vs new Orthogroups.tsv (root N0; lumps genes duplicated between the two roots)", lab_root)):
        r = compare_partitions(universe, lab_old, lab)
        res[key] = r
        say("")
        say(title)
        say(f"  ARI {r['ari']:.3f} | pairs kept {100 * r['kept']:.1f}% | pairs new {100 * r['new']:.1f}%")
        say(f"  groups of 2+ genes: old {r['n_groups_a']}, new {r['n_groups_b']}")
        tot = sum(r["fates"].values())
        for f in ("identical", "merged", "split", "reshuffled"):
            say(f"    {f:<11s} {r['fates'][f]:>7d} old groups ({100 * r['fates'][f] / tot if tot else 0:5.1f}%)  {r['fate_genes'][f]:>8d} genes")
        say(f"  unassigned before, in a group now: {r['joined']} genes; in a group before, unassigned now: {r['left']}")

    with open(os.path.join(a.out, "old_group_fates.tsv"), "w") as fh:
        fh.write("old_group\tn_genes\tfate_at_" + node + "\tn_new_groups_at_" + node + "\tlargest_share_at_" + node +
                 "\tmain_new_hog\tmain_new_hog_n_genes\tfate_at_root\tn_new_groups_at_root\tmain_new_og\tmain_new_og_n_genes\n")
        for grp, (n, fate, nt, share, top, top_n) in sorted(res["A"]["rows"].items()):
            nb, fb, ntb, shb, topb, topb_n = res["B"]["rows"][grp]
            top = top if isinstance(top, str) else "-"
            fh.write(f"{grp}\t{n}\t{fate}\t{nt}\t{share}\t{top}\t{top_n}\t{fb}\t{ntb}\t{topb}\t{topb_n}\n")

    # ---- C orthologue pairs
    say("")
    say("C  orthologue pairs between the shared species (from the gene trees)")
    tot_old = tot_new = tot_both = 0
    rows = []
    for s1, s2 in itertools.combinations(shared, 2):
        po, pn = orthologue_pairs(a.old, s1, s2), orthologue_pairs(a.new, s1, s2)
        if po is None or pn is None:
            continue
        po = {p for p in po if p[0] in universe_set and p[1] in universe_set}
        pn = {p for p in pn if p[0] in universe_set and p[1] in universe_set}
        both = len(po & pn)
        tot_old, tot_new, tot_both = tot_old + len(po), tot_new + len(pn), tot_both + both
        union = len(po) + len(pn) - both
        rows.append((s1, s2, len(po), len(pn), both, both / union if union else 1.0))
    if rows:
        union = tot_old + tot_new - tot_both
        say(f"  {len(rows)} species pairs: old {tot_old} pairs, new {tot_new}; kept {100 * tot_both / tot_old:.1f}% of the old, "
            f"{100 * (tot_new - tot_both) / tot_new:.1f}% of the new ones are new; Jaccard {tot_both / union:.3f}")
        lo = sorted(rows, key=lambda r: r[5])[:3]
        say("  least similar species pairs: " + "; ".join(f"{r[0]} x {r[1]} {r[5]:.3f}" for r in lo))
        with open(os.path.join(a.out, "orthologue_pairs.tsv"), "w") as fh:
            fh.write("species_1\tspecies_2\told_pairs\tnew_pairs\tboth\tjaccard\n")
            for r in rows:
                fh.write("\t".join(map(str, r[:5])) + f"\t{r[5]:.4f}\n")
    else:
        say("  (no pairwise orthologue files found in one of the runs)")

    # ---- named genes
    if a.gene:
        say("")
        say("Named genes: per-species gene counts of their group (shared species; then the added species at the root)")
        added = [s for s in new_sp if s not in set(shared)]
        with open(os.path.join(a.out, "gene_lookups.tsv"), "w") as fh:
            fh.write("name\tquery\tgene\tlevel\tgroup\t" + "\t".join(new_sp) + "\n")
            for spec in a.gene:
                name, _, rest = spec.partition("=")
                sp, _, q = rest.partition(":")
                match = gene_matcher(q)
                hits = sorted(g for g in old_by_sp.get(sp, ()) if match(g))
                if not hits:
                    say(f"  {name} ({sp}:{q}): not found in the old run")
                    continue
                for g in hits[:3]:
                    levels = [("old", old_groups[old_g2[g]], old_g2[g])]
                    if g in new_g2:
                        if not isinstance(lab_node.get(g), tuple):
                            hid = lab_node[g]
                            levels.append((f"new {node}", node_groups[hid], hid))
                        else:
                            levels.append((f"new {node}", {}, "(in no HOG there)"))
                        levels.append(("new root", new_groups[new_g2[g]], new_g2[g]))
                    say(f"  {name} = {g}")
                    for lev, grp, gid in levels:
                        counts = {s: len(grp.get(s, [])) for s in new_sp}
                        sh = "  ".join(f"{code(s)} {counts[s]}" for s in shared)
                        ad = "  ".join(f"{code(s)} {counts[s]}" for s in added) if lev == "new root" else ""
                        say(f"    {lev:<9s} {gid}: {sh}" + (f"  | added: {ad}" if ad else ""))
                        fh.write(f"{name}\t{sp}:{q}\t{g}\t{lev}\t{gid}\t" + "\t".join(str(counts[s]) for s in new_sp) + "\n")
                if len(hits) > 3:
                    say(f"    ({len(hits) - 3} more genes match {q})")

    with open(os.path.join(a.out, "summary.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    say("")
    say(f"written: {a.out}/summary.txt, old_group_fates.tsv, orthologue_pairs.tsv" + (", gene_lookups.tsv" if a.gene else ""))


if __name__ == "__main__":
    main()
