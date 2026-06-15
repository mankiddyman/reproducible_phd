#!/usr/bin/env python3
"""Longest-isoform-per-gene from an AGAT -p protein FASTA.
Headers: '>transcriptID gene=GENEID seq_id=... type=cds'. Keeps, per gene=GENEID,
the single longest protein. Output headers are the GENE id (OrthoFinder genes ==
genes, not isoforms). Distinct gene ids (incl. haplotype copies h1/h2) all kept."""
import sys, re
from collections import OrderedDict

inp, outp = sys.argv[1], sys.argv[2]
gene_re = re.compile(r"gene=(\S+)")

recs = []
hdr = gene = None
seq = []
def flush():
    if hdr is not None:
        recs.append((gene, "".join(seq)))
with open(inp) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush()
            hdr = line[1:]
            m = gene_re.search(hdr)
            gene = m.group(1) if m else hdr.split()[0]
            seq = []
        else:
            seq.append(line)
flush()

best = OrderedDict()  # gene -> (len, seq)
for gene, s in recs:
    if gene not in best or len(s) > best[gene][0]:
        best[gene] = (len(s), s)

with open(outp, "w") as out:
    for gene, (ln, s) in best.items():
        out.write(f">{gene}\n")
        for i in range(0, len(s), 60):
            out.write(s[i:i+60] + "\n")

sys.stderr.write(f"{inp}: {len(recs)} proteins -> {len(best)} genes\n")
