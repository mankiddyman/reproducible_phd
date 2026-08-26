#!/bin/bash
a=$1
[ -s wgd/tre/$a.tre ] && exit 0
seqkit grep -f wgd/ids/$a.ids wgd/all_pep.fa > wgd/fa/$a.fa 2>/dev/null
n=$(grep -c '^>' wgd/fa/$a.fa)
[ "$n" -lt 4 ] && exit 0
mafft --localpair --maxiterate 1000 --quiet --thread 1 wgd/fa/$a.fa > wgd/aln/$a.aln 2>/dev/null
FastTree -quiet -lg wgd/aln/$a.aln > wgd/tre/$a.tre 2>/dev/null
