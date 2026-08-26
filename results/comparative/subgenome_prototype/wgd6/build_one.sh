#!/bin/bash
a=$1
[ -s wgd6/tre/$a.tre ] && exit 0
seqkit grep -f wgd6/ids/$a.ids wgd/all_pep.fa > wgd6/fa/$a.fa
n=$(grep -c '^>' wgd6/fa/$a.fa)
[ "$n" -lt 4 ] && { echo "SKIP $a: only $n seqs" >&2; exit 0; }
mafft --localpair --maxiterate 1000 --quiet --thread 1 wgd6/fa/$a.fa > wgd6/aln/$a.aln
FastTree -quiet -lg wgd6/aln/$a.aln > wgd6/tre/$a.tre
