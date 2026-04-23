#!/usr/bin/env python3
"""
Compile a per-species decisions.yaml from upstream QC outputs.

Stage 2 scope: inputs are the CSV row and the GenomeScope2 summary.txt.
Later stages will extend this script to also consume BlobToolKit, per-role
assembly stats, compleasm summaries, etc., adding fields to the yaml.

Schema rules:
- Numeric estimates with GenomeScope min/max bounds are stored as
  {"min": <int|float>, "max": <int|float>}.
- Derived booleans / strategy labels are top-level scalars.
- Placeholders for fields not yet populated by this stage are written as
  null, so downstream stages can fill them in without schema surprises.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import yaml

# -------- Thresholds / tunables (Stage 2) --------
HETEROZYGOSITY_THRESHOLD_PCT = 1.5  # above this => is_heterozygous: true


def parse_genomescope_summary(path: Path) -> dict:
    """Parse a GenomeScope2 summary.txt into a structured dict.

    The summary has a header section (GenomeScope version, input, output,
    p, k) followed by a table with columns: property, min, max.
    Numeric values may contain commas (e.g. '1,398,083,161 bp') and
    percentages (e.g. '2.31086%'). Strip both before casting.
    """
    text = path.read_text()

    # Header fields (p and k appear as `p = 2` / `k = 21`)
    def find_scalar(key: str) -> str | None:
        m = re.search(rf"^{re.escape(key)}\s*=\s*(\S+)", text, re.MULTILINE)
        return m.group(1) if m else None

    ploidy_str = find_scalar("p")
    k_str = find_scalar("k")

    # Table rows. Each row: <property words>  <min>  <max>
    # The property name can have spaces; min/max are the last two tokens.
    def find_row(property_prefix: str) -> tuple[str, str] | None:
        for line in text.splitlines():
            if line.startswith(property_prefix):
                parts = line.split()
                # Last two whitespace-separated columns are min, max
                # (but values with spaces like "1,398,083,161 bp" confuse this).
                # The table uses multi-space separation; split on 2+ spaces.
                cols = re.split(r"\s{2,}", line.strip())
                if len(cols) >= 3:
                    return cols[-2].strip(), cols[-1].strip()
        return None

    def clean_percent(s: str) -> float:
        return float(s.replace("%", "").replace(",", "").strip())

    def clean_bp(s: str) -> int:
        # "1,398,083,161 bp" -> 1398083161
        return int(s.replace(",", "").replace("bp", "").strip())

    result: dict = {
        "ploidy_from_genomescope": int(ploidy_str) if ploidy_str else None,
        "k": int(k_str) if k_str else None,
    }

    het_row = find_row("Heterozygous")
    if het_row:
        result["heterozygosity_pct"] = {
            "min": clean_percent(het_row[0]),
            "max": clean_percent(het_row[1]),
        }

    homo_row = find_row("Homozygous")
    if homo_row:
        result["homozygosity_pct"] = {
            "min": clean_percent(homo_row[0]),
            "max": clean_percent(homo_row[1]),
        }

    size_row = find_row("Genome Haploid Length")
    if size_row:
        result["genome_haploid_length_bp"] = {
            "min": clean_bp(size_row[0]),
            "max": clean_bp(size_row[1]),
        }

    repeat_row = find_row("Genome Repeat Length")
    if repeat_row:
        result["genome_repeat_length_bp"] = {
            "min": clean_bp(repeat_row[0]),
            "max": clean_bp(repeat_row[1]),
        }

    unique_row = find_row("Genome Unique Length")
    if unique_row:
        result["genome_unique_length_bp"] = {
            "min": clean_bp(unique_row[0]),
            "max": clean_bp(unique_row[1]),
        }

    fit_row = find_row("Model Fit")
    if fit_row:
        result["model_fit_pct"] = {
            "min": clean_percent(fit_row[0]),
            "max": clean_percent(fit_row[1]),
        }

    err_row = find_row("Read Error Rate")
    if err_row:
        result["read_error_rate_pct"] = {
            "min": clean_percent(err_row[0]),
            "max": clean_percent(err_row[1]),
        }

    return result


def derive_decisions(
    species: str,
    exp_ploidy: int,
    centromere: str,
    chr_number_2n: int,
    notes: str,
    gs: dict,
) -> dict:
    """Apply Stage 2 decision logic.

    Current rules:
    - is_heterozygous: True iff heterozygosity_pct.max > HETEROZYGOSITY_THRESHOLD_PCT
    - assembly_strategy:
        - exp_ploidy == 2 -> "phased_diploid"
        - exp_ploidy >  2 -> "phased_polyploid_ragtag"
        - exp_ploidy <  2 -> "primary_only"  (haploid / unusual case)
    - All later-stage fields are null placeholders.
    """
    het_max = None
    if gs.get("heterozygosity_pct"):
        het_max = gs["heterozygosity_pct"]["max"]

    is_het = het_max is not None and het_max > HETEROZYGOSITY_THRESHOLD_PCT

    if exp_ploidy == 2:
        strategy = "phased_diploid"
    elif exp_ploidy > 2:
        strategy = "phased_polyploid_ragtag"
    else:
        strategy = "primary_only"

    decisions = {
        "species": species,
        "csv_inputs": {
            "exp_ploidy": int(exp_ploidy),
            "centromere": centromere,
            "chr_number_2n": int(chr_number_2n),
            "notes": notes or None,
        },
        "genome_estimates": gs,
        "derived": {
            "is_heterozygous": bool(is_het),
            "heterozygosity_threshold_pct": HETEROZYGOSITY_THRESHOLD_PCT,
            "assembly_strategy": strategy,
        },
        # Placeholders for later stages — downstream rules will fill these.
        "contamination": {
            "initial_assessment": None,
            "contaminant_contigs": None,
            "level": None,
        },
        "recommended_assembly": None,
    }
    return decisions


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--exp-ploidy", type=int, required=True)
    ap.add_argument("--centromere", default="")
    ap.add_argument("--chr-number-2n", type=int, required=True)
    ap.add_argument("--notes", default="")
    ap.add_argument("--genomescope-summary", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    if not args.genomescope_summary.is_file():
        print(
            f"ERROR: genomescope summary not found: {args.genomescope_summary}",
            file=sys.stderr,
        )
        return 1

    gs = parse_genomescope_summary(args.genomescope_summary)
    decisions = derive_decisions(
        species=args.species,
        exp_ploidy=args.exp_ploidy,
        centromere=args.centromere,
        chr_number_2n=args.chr_number_2n,
        notes=args.notes,
        gs=gs,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as f:
        yaml.safe_dump(decisions, f, sort_keys=False, default_flow_style=False)

    return 0


if __name__ == "__main__":
    sys.exit(main())
