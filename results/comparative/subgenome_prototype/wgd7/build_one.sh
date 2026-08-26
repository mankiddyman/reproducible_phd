#!/bin/bash
a=$1
[ -s wgd7/tre/$a.tre ] && exit 0
seqkit grep -f wgd7/ids/$a.ids wgd/all_pep.fa > wgd7/fa/$a.fa
n=$(grep -c '^>' wgd7/fa/$a.fa)
[ "$n" -lt 4 ] && { echo "SKIP $a: only $n seqs" >&2; exit 0; }
mafft --localpair --maxiterate 1000 --quiet --thread 1 wgd7/fa/$a.fa > wgd7/aln/$a.aln
FastTree -quiet -lg wgd7/aln/$a.aln > wgd7/tre/$a.tre
