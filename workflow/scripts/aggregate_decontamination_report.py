#!/usr/bin/env python3
"""Aggregate decontamination audit yamls + before/after compleasm into a
single HTML report + TSV summary."""
import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path
import yaml

REPO = Path(".")


def fmt_bp(n):
    if n is None: return "—"
    if n >= 1e9: return f"{n/1e9:.2f}Gb"
    if n >= 1e6: return f"{n/1e6:.1f}Mb"
    if n >= 1e3: return f"{n/1e3:.1f}kb"
    return f"{n}bp"


def parse_compleasm(path):
    """Return dict of S, D, F, M percentages from compleasm summary."""
    if not Path(path).is_file():
        return {}
    txt = Path(path).read_text()
    result = {}
    for key, label in [("S", "Single"), ("D", "Duplicated"), ("F", "Fragmented"), ("M", "Missing")]:
        m = re.search(rf"{key}:(\d+\.\d+)%", txt)
        if m:
            result[label] = float(m.group(1))
    return result


def find_all_audits():
    audits = []
    for path in REPO.glob("results/*/assembly/decontaminated/*/decontamination_audit.yaml"):
        with path.open() as f:
            audit = yaml.safe_load(f)
        sp = audit["species"]
        ass = audit["assembly"]
        # Compleasm before/after
        bf = REPO / f"results/{sp}/qc/compleasm/initial/{ass}/summary.txt"
        af = REPO / f"results/{sp}/qc/compleasm/decontaminated/{ass}/summary.txt"
        # Initial only exists for hap1/hap2 in the current QC scheme;
        # for polyploids targeting p_utg, there's no initial QC at p_utg level.
        audit["compleasm_initial"] = parse_compleasm(bf)
        audit["compleasm_decontaminated"] = parse_compleasm(af)
        # Top dropped classes
        drop_tsv = REPO / f"results/{sp}/assembly/decontaminated/{ass}/drop_contigs.tsv"
        top_dropped = defaultdict(int)
        if drop_tsv.is_file():
            with drop_tsv.open() as f:
                next(f)  # header
                for line in f:
                    parts = line.rstrip().split("\t")
                    if len(parts) >= 3:
                        top_dropped[parts[2]] += int(parts[1])
        audit["top_dropped_classes"] = dict(sorted(top_dropped.items(), key=lambda kv: -kv[1])[:8])
        audits.append(audit)
    audits.sort(key=lambda a: (a["species"], a["assembly"]))
    return audits


def write_tsv(audits, path):
    fields = [
        "species", "assembly",
        "total_contigs", "kept_contigs", "dropped_contigs",
        "total_span", "kept_span", "dropped_span",
        "kept_span_pct",
        "compleasm_initial_S", "compleasm_initial_D",
        "compleasm_decontaminated_S", "compleasm_decontaminated_D",
        "effective_keep_classes",
        "config_hash", "git_commit",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("\t".join(fields) + "\n")
        for a in audits:
            r = a["results"]
            ci = a["compleasm_initial"]
            cd = a["compleasm_decontaminated"]
            row = [
                a["species"], a["assembly"],
                r["total_contigs"], r["kept_contigs"], r["dropped_contigs"],
                r["total_span"], r["kept_span"], r["dropped_span"],
                f"{r['kept_span_pct']:.2f}",
                ci.get("Single", ""), ci.get("Duplicated", ""),
                cd.get("Single", ""), cd.get("Duplicated", ""),
                ",".join(a["filter"]["effective_keep_classes"]),
                a["config_hash"], a["git_commit"],
            ]
            f.write("\t".join(str(x) for x in row) + "\n")


def write_html(audits, path):
    rows = ["""<!DOCTYPE html><html><head><meta charset="utf-8">
    <title>Decontamination report</title>
    <style>
      body { font-family: -apple-system, sans-serif; max-width: 1400px; margin: 2em auto; padding: 0 1em; }
      h1 { border-bottom: 2px solid #333; }
      h2 { background: #f0f0f0; padding: 0.5em; margin-top: 2em; }
      table { border-collapse: collapse; width: 100%; margin: 0.5em 0 1.5em; }
      th, td { padding: 0.3em 0.6em; border-bottom: 1px solid #ddd; text-align: left; font-size: 0.9em; }
      th { background: #eee; }
      .num { text-align: right; font-variant-numeric: tabular-nums; }
      .compleasm-improved { color: #2a7; font-weight: bold; }
      .compleasm-degraded { color: #d33; font-weight: bold; }
      .filter-box { background: #ffd; padding: 0.5em 1em; border: 1px solid #cc9; margin: 0.5em 0; }
    </style></head><body>"""]
    rows.append("<h1>Decontamination report</h1>")
    rows.append(f"<p>{len(audits)} (species, assembly) combinations processed.</p>")

    # Overview table
    rows.append("<h2>Overview</h2>")
    rows.append("<table><thead><tr>"
                "<th>species</th><th>assembly</th>"
                "<th class='num'>contigs (kept/total)</th>"
                "<th class='num'>span (kept/total)</th>"
                "<th class='num'>kept span %</th>"
                "<th class='num'>BUSCO before (S/D)</th>"
                "<th class='num'>BUSCO after (S/D)</th>"
                "<th>filter</th>"
                "</tr></thead><tbody>")
    for a in audits:
        r = a["results"]
        ci = a["compleasm_initial"]
        cd = a["compleasm_decontaminated"]
        bef = f"{ci.get('Single', '—')}/{ci.get('Duplicated', '—')}" if ci else "—"
        aft = f"{cd.get('Single', '—')}/{cd.get('Duplicated', '—')}" if cd else "—"
        # Highlight if BUSCO went down
        cls = ""
        if ci and cd:
            delta = (cd.get('Single', 0) + cd.get('Duplicated', 0)) - (ci.get('Single', 0) + ci.get('Duplicated', 0))
            cls = "compleasm-improved" if delta >= 0 else "compleasm-degraded"
        filt = ", ".join(a["filter"]["effective_keep_classes"])
        rows.append(
            f"<tr><td>{a['species']}</td><td>{a['assembly']}</td>"
            f"<td class='num'>{r['kept_contigs']:,}/{r['total_contigs']:,}</td>"
            f"<td class='num'>{fmt_bp(r['kept_span'])}/{fmt_bp(r['total_span'])}</td>"
            f"<td class='num'>{r['kept_span_pct']:.1f}</td>"
            f"<td class='num'>{bef}</td>"
            f"<td class='num {cls}'>{aft}</td>"
            f"<td>{filt}</td></tr>"
        )
    rows.append("</tbody></table>")

    # Per-species details
    for a in audits:
        rows.append(f"<h2>{a['species']} / {a['assembly']}</h2>")
        rows.append("<div class='filter-box'>")
        rows.append(f"Filter: <strong>{', '.join(a['filter']['effective_keep_classes'])}</strong>")
        if a["filter"]["force_drop_contigs"]:
            rows.append(f"<br>Force-dropped: {len(a['filter']['force_drop_contigs'])} contigs")
        if a["filter"]["force_keep_contigs"]:
            rows.append(f"<br>Force-kept: {len(a['filter']['force_keep_contigs'])} contigs")
        rows.append(f"<br>config_hash: <code>{a['config_hash']}</code> git_commit: <code>{a['git_commit'][:10]}</code>")
        rows.append("</div>")

        # Top dropped classes
        if a["top_dropped_classes"]:
            rows.append("<table><thead><tr><th>top dropped class</th><th class='num'>span</th></tr></thead><tbody>")
            for ccls, span in a["top_dropped_classes"].items():
                rows.append(f"<tr><td>{ccls}</td><td class='num'>{fmt_bp(span)}</td></tr>")
            rows.append("</tbody></table>")

    rows.append("</body></html>")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(rows))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--html", required=True)
    ap.add_argument("--tsv", required=True)
    args = ap.parse_args()

    audits = find_all_audits()
    if not audits:
        print("No audit files found.", file=sys.stderr)
        return 1
    write_tsv(audits, args.tsv)
    write_html(audits, args.html)
    print(f"Wrote {args.tsv} and {args.html} from {len(audits)} audits.", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main() or 0)
