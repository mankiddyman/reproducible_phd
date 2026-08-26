#!/bin/bash
a=$1
[ -s full14/tre/$a.treefile ] && exit 0
seqkit grep -f full14/ids/$a.ids wgd/all_pep.fa > full14/fa/$a.fa 2>/dev/null
grep -c '^>' full14/fa/$a.fa | awk '$1<5{exit 1}' || exit 0
mafft --localpair --maxiterate 1000 --quiet --thread 1 full14/fa/$a.fa > full14/aln/$a.aln 2>/dev/null
iqtree2 -s full14/aln/$a.aln -m MFP -B 1000 -T 1 --prefix full14/tre/$a -quiet -redo >/dev/null 2>&1
