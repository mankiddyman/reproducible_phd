#!/usr/bin/env bash
set -euo pipefail
BASE="${SUBG_BASE:?}"; OUT="$BASE/trees"; J="${J:-24}"
IQ="$(command -v iqtree2)"
n=$(ls "$OUT"/aln/*.fna | wc -l)
mapfile -t picks < <(ls -S "$OUT"/aln/*.fna | sed -n "1p;$((n/2))p;${n}p" | xargs -n1 basename | sed 's/\.fna$//')
lab=(largest median smallest); tot=0; i=0
for a in "${picks[@]}"; do
  L=$(awk 'NR>1 && /^>/ {exit} NR>1 {n+=length($0)} END {print n}' "$OUT/aln/$a.fna")
  T=$(grep -c '^>' "$OUT/aln/$a.fna")
  for s in cds p12 aa; do
    case $s in
      cds) args=(-s "$OUT/aln/$a.fna" -p "$OUT/part/$a.cds.txt" -m GTR+F+G4 --prefix "$OUT/tre_cds/$a") ;;
      p12) args=(-s "$OUT/aln/$a.fna" -p "$OUT/part/$a.p12.txt" -m GTR+F+G4 --prefix "$OUT/tre_p12/$a") ;;
      aa)  [ -f "$OUT/aa/$a.faa" ] || continue
           args=(-s "$OUT/aa/$a.faa" -m LG+F+G4 --prefix "$OUT/tre_aa/$a") ;;
    esac
    t0=$(date +%s.%N); "$IQ" "${args[@]}" -B 1000 -T 1 --quiet -redo >/dev/null 2>&1 \
      || { printf "  %-8s %-4s FAILED\n" "${lab[$i]}" "$s"; continue; }
    t1=$(date +%s.%N); d=$(echo "$t1-$t0" | bc)
    printf "  %-8s %-4s  %5d nt %3d tips  %6.1f core-s\n" "${lab[$i]}" "$s" "$L" "$T" "$d"
    [ "$s" = cds ] && tot=$(echo "$tot+$d" | bc)
  done
  i=$((i+1))
done
mean=$(echo "$tot/3" | bc -l)
printf "\n  CDS mean over the 3: %.1f core-s/tree\n" "$mean"
printf "  projected: %.1f core-h  ->  ~%.0f min wall at -j %s\n" \
  "$(echo "$mean*$n/3600" | bc -l)" "$(echo "$mean*$n/60/$J" | bc -l)" "$J"
printf "  (mean of largest/median/smallest overestimates — size is right-skewed)\n"
