#!/usr/bin/env python3
"""
Normalize a collaborator Helixer GTF so it agrees with its paired FASTA.

Both fixes are AUTO-DETECTED and are no-ops when unnecessary (collaborator
output is inconsistent between haplotypes: Spondias hap1 GTF uses 'scaffold_N'
seqnames while hap2 already uses 'N').

  1. Seqname mismatch: if GTF seqnames are absent from the FASTA but stripping a
     'scaffold_' prefix makes them match, strip it. Aborts if still unmatched.
  2. transcript_id: Helixer wrote a functional DESCRIPTION into transcript_id
     (non-unique across genes). The real transcript identifier is:
       - mRNA lines:         the ID attribute
       - exon/CDS/UTR lines: the Parent attribute (their own ID is feature-level)

gene_id is collaborator-native and is NEVER touched, so downstream functional
annotation joins straight back to their gene names.

Usage: normalize_helixer_gtf.py --gtf IN.gtf --fasta GENOME.fa --out OUT.gtf
"""
import argparse, re, sys
from pathlib import Path

def get_attr(attr, key):
    m = re.search(rf'{key} "([^"]*)"', attr)
    return m.group(1) if m else None

def fasta_names(fasta):
    fai = Path(str(fasta) + ".fai")
    if fai.is_file():
        return {l.split("\t")[0] for l in fai.read_text().splitlines() if l.strip()}
    names = set()
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                names.add(line[1:].split()[0])
    return names

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    fa = fasta_names(a.fasta)
    rows = []
    seq_seen = set()
    with open(a.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                rows.append(line); continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9:
                rows.append(line); continue
            seq_seen.add(p[0]); rows.append(p)

    # --- decide seqname policy ---
    strip = False
    if seq_seen <= fa:
        sys.stderr.write("seqnames: already match FASTA, no change\n")
    else:
        stripped = {s[len("scaffold_"):] if s.startswith("scaffold_") else s
                    for s in seq_seen}
        if stripped <= fa:
            strip = True
            sys.stderr.write("seqnames: stripping 'scaffold_' prefix to match FASTA\n")
        else:
            sys.exit(f"ABORT: GTF seqnames do not match FASTA and stripping "
                     f"'scaffold_' does not fix it.\n"
                     f"  GTF sample:   {sorted(seq_seen)[:5]}\n"
                     f"  FASTA sample: {sorted(fa)[:5]}")

    n_seq = n_tid = 0
    with open(a.out, "w") as out:
        for r in rows:
            if isinstance(r, str):
                out.write(r); continue
            if strip and r[0].startswith("scaffold_"):
                r[0] = r[0][len("scaffold_"):]; n_seq += 1
            attr = r[8]
            if get_attr(attr, "transcript_id") is not None:
                new = get_attr(attr, "ID") if r[2] == "mRNA" else get_attr(attr, "Parent")
                if new is None:
                    g = get_attr(attr, "gene_id")
                    new = f"{g}.1" if g else None
                if new:
                    attr = re.sub(r'transcript_id "[^"]*"', f'transcript_id "{new}"', attr)
                    n_tid += 1
            r[8] = attr
            out.write("\t".join(r) + "\n")

    sys.stderr.write(f"seqnames renamed: {n_seq}\ntranscript_id set: {n_tid}\n")

if __name__ == "__main__":
    main()
