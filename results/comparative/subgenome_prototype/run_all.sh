#!/usr/bin/env bash
# Dionaea homeolog divergence relative to Nepenthes, genome-wide.
#   stage 1  triplets      1 Nep(dom) : 2 Dio in different tandem arrays
#   stage 2  align         mafft L-INS-i on peptides -> pal2nal back-translation
#   stage 3  distances     protein p-distance + Nei-Gojobori dN/dS, all 3 pairs
#   stage 4  plots         R
# usage: bash run_all.sh [outdir]
set -euo pipefail
OUT="${1:-work}"
RP=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
G=$RP/results/comparative/genespace/results
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
mkdir -p "$OUT/aln"
module load seqkit 2>/dev/null || true

# ---------- stage 1 ----------
if [ ! -s triplets_all.tsv ]; then
  echo "[1] building triplets"
  awk -F'\t' 'NR>1 && ($8=="Nepenthes_gracilis" || $8=="Dionaea_muscipula") {
      print $12"\t"$8"\t"$1"\t"$4"\t"$9 }' $G/combBed.txt \
  | awk -F'\t' '
    $2=="Nepenthes_gracilis" && $3 ~ /_dom$/ { nep[$1]=$4; nn[$1]++ }
    $2=="Dionaea_muscipula" { dn[$1]++; id[$1]=id[$1]" "$4; ar[$1]=ar[$1]" "$5; ch[$1]=ch[$1]" "$3 }
    END { print "hog\tnep\tdio1\tdio2\tchr1\tchr2\tarrangement"
      for (k in nep) if (nn[k]==1 && dn[k]==2) {
        split(id[k],I," "); split(ar[k],A," "); split(ch[k],C," ")
        if (A[1]==A[2]) continue
        print k"\t"nep[k]"\t"I[1]"\t"I[2]"\t"C[1]"\t"C[2]"\t"(C[1]==C[2] ? "same_chr_translocated":"different_chr")
      }}' > triplets_all.tsv
fi
echo "[1] triplets: $(($(wc -l < triplets_all.tsv)-1))"

# ---------- sequences ----------
if [ ! -s all_pep.faa ]; then
  cat $RP/results/comparative/genespace/peptide/Nepenthes_gracilis.fa \
      $RP/results/comparative/genespace/peptide/Dionaea_muscipula.fa > all_pep.faa
fi
if [ ! -s all_cds.fna ]; then
  for sp in Nepenthes_gracilis Dionaea_muscipula; do
    gffread $RP/results/$sp/annotation/final/$sp.final.gff3 \
      -g $RP/results/$sp/assembly_final/external_collapsed/${sp}_chr.fa \
      -x - 2>/dev/null | sed '/^>/! s/\.//g' >> all_cds.fna
  done
fi

# ---------- stage 2 ----------
echo "[2] aligning"
: > "$OUT/aln_index.tsv"
tail -n +2 triplets_all.tsv | while IFS=$'\t' read -r hog nep d1 d2 c1 c2 arr; do
  hid=$(echo "$hog" | tr ' ' '_' | cut -d_ -f1)
  d="$OUT/aln/$hid"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$hid" "$nep" "$d1" "$d2" "$c1" "$c2" "$arr" >> "$OUT/aln_index.tsv"
  [ -s "$d.pal2nal.fna" ] && continue
  printf "%s\n%s\n%s\n" "$nep" "$d1" "$d2" > "$d.ids"
  seqkit grep -f "$d.ids" all_pep.faa 2>/dev/null > "$d.faa"
  seqkit grep -f "$d.ids" all_cds.fna 2>/dev/null > "$d.fna"
  [ "$(grep -c '^>' "$d.faa")" -eq 3 ] || continue
  mafft --maxiterate 1000 --localpair --quiet "$d.faa" > "$d.aln.faa" 2>/dev/null
  pal2nal.pl "$d.aln.faa" "$d.fna" -output fasta -nogap > "$d.pal2nal.fna" 2>/dev/null || true
done
echo "[2] aligned: $(ls "$OUT"/aln/*.pal2nal.fna 2>/dev/null | wc -l)"

echo "[3] distances";  python3 distances.py "$OUT"
echo "[4] plots";      Rscript plot_distances.R
