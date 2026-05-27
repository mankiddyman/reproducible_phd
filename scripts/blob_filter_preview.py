#!/usr/bin/env python3
"""Quick-and-dirty filter preview for blobtoolkit outputs.

For each blobdir, read the per-contig category assignments at multiple
taxonomic ranks and report:
  - Total contigs / span
  - What "keep Magnoliopsida + no-hit" retains and what it drops
  - Top hits at each taxonomic rank among contigs we would drop
  - Top hits at each rank among contigs we would keep

Usage:
  python3 scripts/blob_filter_preview.py
  python3 scripts/blob_filter_preview.py --species Drosera_paradoxa --assembly hap1
  python3 scripts/blob_filter_preview.py --html out.html
"""
import argparse
import gzip
import json
from pathlib import Path
from collections import Counter, defaultdict
import sys

REPO = Path("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd")

# Filter logic: keep contigs whose buscoregions_class is Magnoliopsida or no-hit.
# Tunable later; this is the proposed default.
KEEP_CLASSES = {"Magnoliopsida", "no-hit"}

RANKS = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]


def load_gz_json(path):
    """Load a .json.gz file or its uncompressed sibling."""
    p = Path(path)
    if p.exists():
        with open(p) as f:
            return json.load(f)
    gz = Path(str(p) + ".gz")
    if gz.exists():
        with gzip.open(gz, "rt") as f:
            return json.load(f)
    return None


def load_category_field(blobdir, field_id):
    """Load a category JSON like buscoregions_class.json[.gz]; returns the
    field dict with .values (list of category indices) and .keys (category names)."""
    return load_gz_json(blobdir / f"{field_id}.json")


def load_variable_field(blobdir, field_id):
    """Load a variable JSON like length.json[.gz]; returns {.values: [...]}."""
    return load_gz_json(blobdir / f"{field_id}.json")


def per_contig_categories(blobdir, rank):
    """Return list of category names per contig at the given rank."""
    field = load_category_field(blobdir, f"buscoregions_{rank}")
    if field is None:
        return None
    keys = field.get("keys", [])
    values = field.get("values", [])
    # keys[i] = name string; values[contig_idx] = key index. Use "no-hit" for missing
    return [keys[v] if 0 <= v < len(keys) else "no-hit" for v in values]


def analyse(blobdir):
    """Return a dict of stats and per-contig classifications for one blobdir."""
    meta = load_gz_json(blobdir / "meta.json")
    if meta is None:
        return None

    identifiers = load_variable_field(blobdir, "identifiers")
    lengths = load_variable_field(blobdir, "length")
    gc = load_variable_field(blobdir, "gc")
    if not (identifiers and lengths):
        return None

    contig_ids = identifiers.get("values", [])
    contig_lens = lengths.get("values", [])
    contig_gc = gc.get("values", []) if gc else [None] * len(contig_ids)

    # Per-rank classification
    classifications = {}
    for rank in RANKS:
        classifications[rank] = per_contig_categories(blobdir, rank)

    # Apply the filter using buscoregions_class
    class_calls = classifications.get("class") or ["no-hit"] * len(contig_ids)
    keep_mask = [c in KEEP_CLASSES for c in class_calls]

    total_n = len(contig_ids)
    total_span = sum(contig_lens)
    keep_n = sum(keep_mask)
    keep_span = sum(l for l, k in zip(contig_lens, keep_mask) if k)

    # Count taxon abundance among DROPPED and KEPT contigs at each rank
    by_rank_drop = {}
    by_rank_keep = {}
    for rank, calls in classifications.items():
        if calls is None:
            continue
        drop_counter = Counter()
        keep_counter = Counter()
        for call, l, k in zip(calls, contig_lens, keep_mask):
            if k:
                keep_counter[call] += l
            else:
                drop_counter[call] += l
        by_rank_drop[rank] = drop_counter
        by_rank_keep[rank] = keep_counter

    return {
        "meta": meta,
        "n_contigs": total_n,
        "total_span": total_span,
        "keep_n": keep_n,
        "keep_span": keep_span,
        "drop_n": total_n - keep_n,
        "drop_span": total_span - keep_span,
        "by_rank_drop": by_rank_drop,
        "by_rank_keep": by_rank_keep,
        "contig_ids": contig_ids,
        "contig_lens": contig_lens,
        "contig_gc": contig_gc,
        "class_calls": class_calls,
    }


def fmt_bp(n):
    if n >= 1e9: return f"{n/1e9:.2f}Gb"
    if n >= 1e6: return f"{n/1e6:.1f}Mb"
    if n >= 1e3: return f"{n/1e3:.1f}kb"
    return f"{n}bp"


def print_report(species, assembly, stats):
    print(f"\n{'='*70}")
    print(f" {species} / {assembly}")
    print(f"{'='*70}")
    print(f" Total: {stats['n_contigs']:,} contigs, {fmt_bp(stats['total_span'])}")
    print(f" KEEP : {stats['keep_n']:,} contigs ({stats['keep_n']/stats['n_contigs']*100:.1f}%), "
          f"{fmt_bp(stats['keep_span'])} ({stats['keep_span']/stats['total_span']*100:.1f}%)")
    print(f" DROP : {stats['drop_n']:,} contigs ({stats['drop_n']/stats['n_contigs']*100:.1f}%), "
          f"{fmt_bp(stats['drop_span'])} ({stats['drop_span']/stats['total_span']*100:.1f}%)")
    print()
    print(f" Top categories among DROPPED contigs (by span):")
    for rank in RANKS:
        if rank not in stats["by_rank_drop"]:
            continue
        top = stats["by_rank_drop"][rank].most_common(5)
        if not top:
            continue
        print(f"   {rank:14s}: " + ", ".join(f"{name} ({fmt_bp(span)})" for name, span in top))
    print()
    print(f" Top categories among KEPT contigs (by span):")
    for rank in RANKS:
        if rank not in stats["by_rank_keep"]:
            continue
        top = stats["by_rank_keep"][rank].most_common(5)
        if not top:
            continue
        print(f"   {rank:14s}: " + ", ".join(f"{name} ({fmt_bp(span)})" for name, span in top))


def write_html(reports, out_path):
    """Write a single HTML page summarising all (species, assembly) reports."""
    rows = []
    rows.append("""<!DOCTYPE html><html><head><meta charset="utf-8">
    <title>Blob filter preview</title>
    <style>
      body { font-family: -apple-system, sans-serif; max-width: 1400px; margin: 2em auto; padding: 0 1em; }
      h1 { border-bottom: 2px solid #333; }
      h2 { background: #f0f0f0; padding: 0.5em; margin-top: 2em; }
      .summary { display: grid; grid-template-columns: repeat(4, 1fr); gap: 1em; margin: 1em 0; }
      .summary div { background: #f8f8f8; padding: 0.8em; border-left: 4px solid #2a7; }
      .summary .drop { border-left-color: #d33; }
      table { border-collapse: collapse; width: 100%; margin: 0.5em 0 1.5em; }
      th, td { padding: 0.3em 0.6em; border-bottom: 1px solid #ddd; text-align: left; font-size: 0.9em; }
      th { background: #eee; }
      .drop-row td:last-child { color: #d33; }
      .keep-row td:last-child { color: #2a7; }
      .num { text-align: right; font-variant-numeric: tabular-nums; }
      .filter { background: #ffd; padding: 0.5em 1em; border: 1px solid #cc9; }
    </style></head><body>""")
    rows.append("<h1>Blobtoolkit filter preview</h1>")
    rows.append(f'<div class="filter">Filter rule: <strong>keep contigs where class ∈ {sorted(KEEP_CLASSES)}</strong></div>')

    # Build cross-species summary table
    rows.append("<h2>Cross-species summary</h2>")
    rows.append("<table><thead><tr><th>species</th><th>assembly</th>"
                "<th class='num'>n contigs</th><th class='num'>total span</th>"
                "<th class='num'>keep n</th><th class='num'>keep %</th>"
                "<th class='num'>keep span</th><th class='num'>keep span %</th>"
                "<th class='num'>drop n</th><th class='num'>drop span</th>"
                "</tr></thead><tbody>")
    for (sp, ass, stats) in reports:
        rows.append(
            f"<tr><td>{sp}</td><td>{ass}</td>"
            f"<td class='num'>{stats['n_contigs']:,}</td>"
            f"<td class='num'>{fmt_bp(stats['total_span'])}</td>"
            f"<td class='num'>{stats['keep_n']:,}</td>"
            f"<td class='num'>{stats['keep_n']/stats['n_contigs']*100:.1f}</td>"
            f"<td class='num'>{fmt_bp(stats['keep_span'])}</td>"
            f"<td class='num'>{stats['keep_span']/stats['total_span']*100:.1f}</td>"
            f"<td class='num'>{stats['drop_n']:,}</td>"
            f"<td class='num'>{fmt_bp(stats['drop_span'])}</td></tr>"
        )
    rows.append("</tbody></table>")

    # Cross-species: aggregate "what is being dropped" at each rank
    rows.append("<h2>Cross-species: most common DROPPED taxa</h2>")
    rows.append("<p>For each taxonomic rank, how many species had this taxon dropped, and total span across all species.</p>")
    for rank in RANKS:
        # taxon -> (n_species_seen, total_span)
        agg = defaultdict(lambda: [0, 0])
        for (sp, ass, stats) in reports:
            seen_here = set()
            for taxon, span in stats["by_rank_drop"].get(rank, {}).items():
                if taxon and taxon != "no-hit":
                    agg[taxon][1] += span
                    if taxon not in seen_here:
                        agg[taxon][0] += 1
                        seen_here.add(taxon)
        if not agg:
            continue
        rows.append(f"<h3>{rank}</h3>")
        rows.append("<table><thead><tr><th>taxon</th><th class='num'>species count</th><th class='num'>total dropped span</th></tr></thead><tbody>")
        sorted_taxa = sorted(agg.items(), key=lambda kv: -kv[1][1])[:15]
        for taxon, (n_sp, span) in sorted_taxa:
            rows.append(f"<tr class='drop-row'><td>{taxon}</td><td class='num'>{n_sp}</td><td class='num'>{fmt_bp(span)}</td></tr>")
        rows.append("</tbody></table>")

    # Per-(species, assembly) detail
    for (sp, ass, stats) in reports:
        rows.append(f"<h2>{sp} / {ass}</h2>")
        rows.append(f"<div class='summary'>")
        rows.append(f"<div><strong>{stats['n_contigs']:,}</strong><br>total contigs</div>")
        rows.append(f"<div><strong>{fmt_bp(stats['total_span'])}</strong><br>total span</div>")
        rows.append(f"<div><strong>{stats['keep_n']:,} ({stats['keep_n']/stats['n_contigs']*100:.1f}%)</strong><br>keep contigs</div>")
        rows.append(f"<div class='drop'><strong>{fmt_bp(stats['drop_span'])} ({stats['drop_span']/stats['total_span']*100:.1f}%)</strong><br>drop span</div>")
        rows.append("</div>")

        rows.append("<table><thead><tr><th>rank</th><th>top DROPPED (by span)</th><th>top KEPT (by span)</th></tr></thead><tbody>")
        for rank in RANKS:
            drop_top = stats["by_rank_drop"].get(rank, Counter()).most_common(5)
            keep_top = stats["by_rank_keep"].get(rank, Counter()).most_common(5)
            drop_txt = "; ".join(f"{n} ({fmt_bp(s)})" for n, s in drop_top) or "—"
            keep_txt = "; ".join(f"{n} ({fmt_bp(s)})" for n, s in keep_top) or "—"
            rows.append(f"<tr><td>{rank}</td><td>{drop_txt}</td><td>{keep_txt}</td></tr>")
        rows.append("</tbody></table>")

    rows.append("</body></html>")
    Path(out_path).write_text("\n".join(rows))
    print(f"\nHTML report written to: {out_path}", file=sys.stderr)


def find_blobdirs(species_filter=None, assembly_filter=None):
    """Walk results/ to find all blobdirs."""
    results = []
    for sp_dir in (REPO / "results").glob("Drosera_*"):
        sp = sp_dir.name
        if species_filter and species_filter not in sp:
            continue
        for ass_dir in (sp_dir / "blobtoolkit/initial").glob("*"):
            assembly = ass_dir.name
            if assembly_filter and assembly_filter != assembly:
                continue
            # blobdir is at: <ass_dir>/output/blobtoolkit/<species>/
            blobdir = ass_dir / "output" / "blobtoolkit" / sp
            if not blobdir.is_dir():
                continue
            if not ((blobdir / "meta.json").exists() or (blobdir / "meta.json.gz").exists()):
                continue
            results.append((sp, assembly, blobdir))
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", help="Filter to species name substring")
    ap.add_argument("--assembly", help="Filter to assembly (hap1, hap2, p_utg)")
    ap.add_argument("--html", help="Write HTML report to this path")
    args = ap.parse_args()

    blobdirs = find_blobdirs(args.species, args.assembly)
    if not blobdirs:
        print("No blobdirs found.")
        return 1

    reports = []
    for sp, ass, blobdir in blobdirs:
        print(f"Analysing {sp} / {ass} ...", file=sys.stderr)
        stats = analyse(blobdir)
        if stats is None:
            print(f"  failed", file=sys.stderr)
            continue
        reports.append((sp, ass, stats))
        if not args.html:
            print_report(sp, ass, stats)

    if args.html:
        write_html(reports, args.html)


if __name__ == "__main__":
    sys.exit(main() or 0)
