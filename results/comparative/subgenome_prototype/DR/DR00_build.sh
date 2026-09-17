#!/bin/bash
# one locus: peptides -> mafft -> pal2nal -> codon alignment
set -u
a=$1
[ -s DR/codon/$a.fna ] && exit 0
seqkit grep -f DR/ids/$a.ids wgd/all_pep.fa        > DR/fa/$a.fa    2>/dev/null
seqkit grep -f DR/ids/$a.ids cds/all_cds_tagged.fa > DR/cdsfa/$a.fna 2>/dev/null
n=$(grep -c '^>' DR/fa/$a.fa 2>/dev/null || echo 0)
[ "$n" -lt 3 ] && exit 0
mafft --localpair --maxiterate 1000 --quiet --thread 1 DR/fa/$a.fa > DR/aln/$a.aln 2>/dev/null
/usr/bin/pal2nal.pl DR/aln/$a.aln DR/cdsfa/$a.fna -output fasta 2>/dev/null > DR/codon/$a.fna
[ -s DR/codon/$a.fna ] || rm -f DR/codon/$a.fna
