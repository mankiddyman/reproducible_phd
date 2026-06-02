#!/usr/bin/env python3
"""Split + rename a frozen assembly using a per-species naming map.

Sequences listed in the map (old_name,new_name) -> renamed, written to
--out-chr in map order (so chr1_hap1..chrN_hapM come out ordered).
Sequences NOT in the map -> renamed unplaced_NNN by descending length,
written to --out-unplaced.
"""
import argparse, csv, sys

def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)

def write_fasta(path, records, width=60):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i+width] + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--out-chr", required=True)
    ap.add_argument("--out-unplaced", required=True)
    args = ap.parse_args()

    # naming map (preserve order for chr output)
    name_map, order = {}, []
    with open(args.map) as fh:
        r = csv.DictReader(fh)
        for row in r:
            name_map[row["old_name"]] = row["new_name"]
            order.append(row["old_name"])

    seqs = dict(read_fasta(args.fasta))

    # chromosomes in map order; error if a mapped scaffold is missing
    missing = [o for o in order if o not in seqs]
    if missing:
        sys.exit(f"ERROR: mapped scaffolds absent from fasta: {missing}")
    chr_records = [(name_map[o], seqs[o]) for o in order]

    # unplaced = everything not in the map, by descending length
    unplaced_names = [n for n in seqs if n not in name_map]
    unplaced_names.sort(key=lambda n: len(seqs[n]), reverse=True)
    unplaced_records = [
        (f"unplaced_{i:03d}", seqs[n]) for i, n in enumerate(unplaced_names, 1)
    ]

    write_fasta(args.out_chr, chr_records)
    write_fasta(args.out_unplaced, unplaced_records)

    print(f"chr scaffolds : {len(chr_records)} -> {args.out_chr}")
    print(f"unplaced      : {len(unplaced_records)} -> {args.out_unplaced}")

if __name__ == "__main__":
    main()
