#!/usr/bin/env bash
# Per-triplet (1 Nepenthes : 2 Dionaea homeologs) pairwise distances.
# usage: bash run_triplets.sh [triplets.tsv] [outdir]
set -euo pipefail
TRIP="${1:-triplets.tsv}"
OUT="${2:-work}"
mkdir -p "$OUT/aln"
module load seqkit 2>/dev/null || true

N=0
: > "$OUT/aln_index.tsv"
tail -n +2 "$TRIP" | while IFS=$'\t' read -r hog nep s5 s9; do
  hid=$(echo "$hog" | tr ' ' '_' | cut -d_ -f1)
  d="$OUT/aln/$hid"
  printf "%s\t%s\t%s\t%s\n" "$hid" "$nep" "$s5" "$s9" >> "$OUT/aln_index.tsv"
  [ -f "$d.pal2nal.fna" ] && continue
  printf "%s\n%s\n%s\n" "$nep" "$s5" "$s9" > "$d.ids"
  seqkit grep -f "$d.ids" all_pep.faa 2>/dev/null > "$d.faa"
  seqkit grep -f "$d.ids" all_cds.fna 2>/dev/null > "$d.fna"
  [ "$(grep -c '^>' "$d.faa")" -eq 3 ] || { echo "SKIP $hid (pep)" >&2; continue; }
  mafft --maxiterate 1000 --localpair --quiet "$d.faa" > "$d.aln.faa" 2>/dev/null
  pal2nal.pl "$d.aln.faa" "$d.fna" -output fasta -nogap > "$d.pal2nal.fna" 2>/dev/null || true
  [ -s "$d.pal2nal.fna" ] || echo "SKIP $hid (pal2nal)" >&2
done
echo "aligned: $(ls "$OUT"/aln/*.pal2nal.fna 2>/dev/null | wc -l) / $(($(wc -l < "$TRIP")-1))"
