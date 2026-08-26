#!/bin/bash
a=$1
[ -s wgd7/codon/$a.fna ] && exit 0
[ -s wgd7/aln/$a.aln ] || exit 0
seqkit grep -f wgd7/ids/$a.ids cds/all_cds_tagged.fa > wgd7/cdsfa/$a.fna 2>/dev/null
np=$(grep -c '^>' wgd7/aln/$a.aln); nc=$(grep -c '^>' wgd7/cdsfa/$a.fna)
[ "$np" -ne "$nc" ] && { echo "SKIP $a: aln=$np cds=$nc" >&2; exit 0; }
/usr/bin/pal2nal.pl wgd7/aln/$a.aln wgd7/cdsfa/$a.fna -output fasta -nogap 2>/dev/null \
  > wgd7/codon/$a.fna
[ -s wgd7/codon/$a.fna ] || { echo "PAL2NAL FAILED $a" >&2; rm -f wgd7/codon/$a.fna; }
