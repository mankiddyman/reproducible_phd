#!/usr/bin/env python3
"""Decontaminate an assembly using blobtoolkit buscoregions_class assignments.

Reads:
  - blobdir directory (with meta.json and *_class.json[.gz] files)
  - decontamination config (yaml)
  - input FASTA (only metadata; sequence extraction happens in the rule via seqkit)

Writes:
  - keep_contigs.tsv (contig_id, length, class, decision_reason)
  - drop_contigs.tsv (same schema)
  - decontamination_audit.yaml (full provenance)
  - stats.tsv (before/after counts and spans)
  - keep_contigs.list (one contig id per line, for seqkit grep)
"""
import argparse
import gzip
import hashlib
import json
import subprocess
import sys
import yaml
from datetime import datetime, timezone
from pathlib import Path


def load_gz_json(path):
    p = Path(path)
    if p.exists():
        with open(p) as f:
            return json.load(f)
    gz = Path(str(p) + ".gz")
    if gz.exists():
        with gzip.open(gz, "rt") as f:
            return json.load(f)
    return None


def load_class_calls(blobdir):
    """Return per-contig class names (using buscoregions_class)."""
    field = load_gz_json(blobdir / "buscoregions_class.json")
    if field is None:
        raise FileNotFoundError(f"buscoregions_class.json[.gz] not found in {blobdir}")
    keys = field.get("keys", [])
    values = field.get("values", [])
    return [keys[v] if 0 <= v < len(keys) else "no-hit" for v in values]


def load_identifiers(blobdir):
    field = load_gz_json(blobdir / "identifiers.json")
    if field is None:
        raise FileNotFoundError(f"identifiers.json[.gz] not found in {blobdir}")
    return field["values"]


def load_lengths(blobdir):
    field = load_gz_json(blobdir / "length.json")
    if field is None:
        raise FileNotFoundError(f"length.json[.gz] not found in {blobdir}")
    return field["values"]


def git_commit():
    """Return current git commit, or 'unknown' if not in a git repo."""
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL
        )
        return out.decode().strip()
    except Exception:
        return "unknown"


def config_hash(cfg_dict):
    """Stable hash of the decontamination config used."""
    payload = json.dumps(cfg_dict, sort_keys=True).encode()
    return hashlib.sha256(payload).hexdigest()[:16]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blobdir", required=True, help="Path to blobdir directory")
    ap.add_argument("--config", required=True, help="config/decontamination.yaml path")
    ap.add_argument("--species", required=True)
    ap.add_argument("--assembly", required=True, help="hap1, hap2, p_utg etc.")
    ap.add_argument("--out-keep-tsv", required=True)
    ap.add_argument("--out-drop-tsv", required=True)
    ap.add_argument("--out-keep-list", required=True, help="One contig id per line, for seqkit")
    ap.add_argument("--out-audit-yaml", required=True)
    ap.add_argument("--out-stats-tsv", required=True)
    args = ap.parse_args()

    blobdir = Path(args.blobdir)
    cfg_path = Path(args.config)

    # Load decontamination config
    if cfg_path.is_file():
        with cfg_path.open() as f:
            cfg = yaml.safe_load(f) or {}
    else:
        cfg = {}
    cfg.setdefault("defaults", {"keep_classes": ["Magnoliopsida", "no-hit"]})
    cfg.setdefault("species_overrides", {})
    if cfg["species_overrides"] is None:
        cfg["species_overrides"] = {}

    defaults = cfg["defaults"]
    overrides = cfg["species_overrides"].get(args.species, {}) or {}

    base_keep = set(defaults.get("keep_classes", []))
    additional_keep = set(overrides.get("additional_keep_classes", []) or [])
    effective_keep_classes = base_keep | additional_keep

    force_drop = set(overrides.get("force_drop_contigs", []) or [])
    force_keep = set(overrides.get("force_keep_contigs", []) or [])

    # Load blobdir data
    contig_ids = load_identifiers(blobdir)
    contig_lens = load_lengths(blobdir)
    class_calls = load_class_calls(blobdir)
    n = len(contig_ids)
    if not (len(contig_lens) == len(class_calls) == n):
        raise RuntimeError(
            f"Length mismatch: ids={n}, lens={len(contig_lens)}, "
            f"classes={len(class_calls)}"
        )

    # Apply filter
    keep_rows = []
    drop_rows = []
    drop_class_span = {}
    for cid, clen, ccls in zip(contig_ids, contig_lens, class_calls):
        if cid in force_keep:
            decision = "keep"
            reason = "force_keep_override"
        elif cid in force_drop:
            decision = "drop"
            reason = "force_drop_override"
        elif ccls in effective_keep_classes:
            decision = "keep"
            reason = f"class={ccls}"
        else:
            decision = "drop"
            reason = f"class={ccls}"

        row = (cid, clen, ccls, reason)
        if decision == "keep":
            keep_rows.append(row)
        else:
            drop_rows.append(row)
            drop_class_span[ccls] = drop_class_span.get(ccls, 0) + clen

    # Write outputs
    Path(args.out_keep_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_keep_tsv, "w") as f:
        f.write("contig_id\tlength\tclass\treason\n")
        for cid, clen, ccls, reason in keep_rows:
            f.write(f"{cid}\t{clen}\t{ccls}\t{reason}\n")

    with open(args.out_drop_tsv, "w") as f:
        f.write("contig_id\tlength\tclass\treason\n")
        for cid, clen, ccls, reason in drop_rows:
            f.write(f"{cid}\t{clen}\t{ccls}\t{reason}\n")

    with open(args.out_keep_list, "w") as f:
        for cid, _, _, _ in keep_rows:
            f.write(f"{cid}\n")

    total_n = n
    keep_n = len(keep_rows)
    drop_n = len(drop_rows)
    total_span = sum(contig_lens)
    keep_span = sum(r[1] for r in keep_rows)
    drop_span = sum(r[1] for r in drop_rows)

    with open(args.out_stats_tsv, "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"total_contigs\t{total_n}\n")
        f.write(f"kept_contigs\t{keep_n}\n")
        f.write(f"dropped_contigs\t{drop_n}\n")
        f.write(f"total_span\t{total_span}\n")
        f.write(f"kept_span\t{keep_span}\n")
        f.write(f"dropped_span\t{drop_span}\n")
        f.write(f"kept_contig_pct\t{keep_n/total_n*100:.2f}\n")
        f.write(f"kept_span_pct\t{keep_span/total_span*100:.2f}\n")

    audit = {
        "species": args.species,
        "assembly": args.assembly,
        "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "git_commit": git_commit(),
        "config_hash": config_hash(cfg),
        "filter": {
            "default_keep_classes": sorted(base_keep),
            "additional_keep_classes": sorted(additional_keep),
            "effective_keep_classes": sorted(effective_keep_classes),
            "force_drop_contigs": sorted(force_drop),
            "force_keep_contigs": sorted(force_keep),
        },
        "results": {
            "total_contigs": total_n,
            "kept_contigs": keep_n,
            "dropped_contigs": drop_n,
            "total_span": total_span,
            "kept_span": keep_span,
            "dropped_span": drop_span,
            "kept_contig_pct": round(keep_n / total_n * 100, 2),
            "kept_span_pct": round(keep_span / total_span * 100, 2),
            "drop_class_span": {k: drop_class_span[k] for k in sorted(drop_class_span)},
        },
    }
    with open(args.out_audit_yaml, "w") as f:
        yaml.safe_dump(audit, f, sort_keys=False)

    print(f"OK: {args.species}/{args.assembly} — kept {keep_n}/{total_n} contigs "
          f"({keep_span/total_span*100:.1f}% of span)", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main() or 0)
