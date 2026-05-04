#!/usr/bin/env python3
"""
Generate a per-species samplesheet for the sanger-tol/blobtoolkit Nextflow pipeline.

Per upstream README: "The first column (sample name) must be unique. If you have
multiple read sets from the same actual sample, make sure you edit the sample
names to make them unique."

We name each read set as {species}_run{N} where N is 1-indexed by file order.
"""

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--species-table", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.species_table).fillna("")
    rows = df[df["species_id"] == args.species]
    if rows.empty:
        print(f"ERROR: species {args.species} not in {args.species_table}", file=sys.stderr)
        return 1

    raw = rows.iloc[0]["genomic_fastq"]
    fastqs = [p.strip() for p in str(raw).split(";") if p.strip()]
    if not fastqs:
        print(f"ERROR: no genomic_fastq paths for {args.species}", file=sys.stderr)
        return 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "datatype", "datafile", "library_layout"])
        for i, fq in enumerate(fastqs, start=1):
            sample_name = f"{args.species}_run{i}"
            writer.writerow([sample_name, "pacbio", fq, "SINGLE"])

    return 0


if __name__ == "__main__":
    sys.exit(main())
