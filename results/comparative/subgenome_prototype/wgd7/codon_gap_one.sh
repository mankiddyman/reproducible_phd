#!/bin/bash
a=$1
[ -s wgd7/codon_gap/$a.fna ] && exit 0
[ -s wgd7/cdsfa/$a.fna ] || exit 0
/usr/bin/pal2nal.pl wgd7/aln/$a.aln wgd7/cdsfa/$a.fna -output fasta 2>/dev/null \
  > wgd7/codon_gap/$a.fna
[ -s wgd7/codon_gap/$a.fna ] || rm -f wgd7/codon_gap/$a.fna
