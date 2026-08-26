#!/bin/bash
a=$1
[ -s wgd7iq/tre/$a.tre ] && exit 0
[ -s wgd7/aln/$a.aln ] || { echo "no aln $a" >&2; exit 0; }
n=$(grep -c '^>' wgd7/aln/$a.aln)
[ "$n" -lt 4 ] && exit 0
iqtree2 -s wgd7/aln/$a.aln -m LG+F+G4 -B 1000 -T 1 \
        --prefix wgd7iq/tre/$a -quiet -redo > /dev/null 2>&1
[ -s wgd7iq/tre/$a.treefile ] && mv wgd7iq/tre/$a.treefile wgd7iq/tre/$a.tre
rm -f wgd7iq/tre/$a.ckp.gz wgd7iq/tre/$a.iqtree wgd7iq/tre/$a.log \
      wgd7iq/tre/$a.mldist wgd7iq/tre/$a.splits.nex wgd7iq/tre/$a.bionj \
      wgd7iq/tre/$a.contree wgd7iq/tre/$a.model.gz wgd7iq/tre/$a.uniqueseq.phy
