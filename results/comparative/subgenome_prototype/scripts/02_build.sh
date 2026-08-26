#!/usr/bin/env bash
# Stage 02 — build.  ./02_build.sh {time|cds|aa|p12}
# One single-threaded IQ-TREE per anchor under GNU parallel. IQ-TREE scales
# badly on ~10-tip / ~500-column alignments, so -T 24 on one small tree wastes
# most of the cores; 24 independent jobs does not.
set -euo pipefail
BASE="${SUBG_BASE:-$PWD}"; OUT="$BASE/trees"; J="${J:-24}"
IQ="$(command -v iqtree2 || command -v iqtree)"; MODE="${1:-time}"

run_one () {
  local set="$1" anchor="$2"
  case "$set" in
    cds) "$IQ" -s "$OUT/aln/$anchor.fna" -p "$OUT/part/$anchor.cds.txt" \
               -m GTR+F+G4 -B 1000 -T 1 --quiet -redo --prefix "$OUT/tre_cds/$anchor" ;;
    p12) "$IQ" -s "$OUT/aln/$anchor.fna" -p "$OUT/part/$anchor.p12.txt" \
               -m GTR+F+G4 -B 1000 -T 1 --quiet -redo --prefix "$OUT/tre_p12/$anchor" ;;
    aa)  "$IQ" -s "$OUT/aa/$anchor.faa" \
               -m LG+F+G4 -B 1000 -T 1 --quiet -redo --prefix "$OUT/tre_aa/$anchor" ;;
  esac
}
export -f run_one; export OUT IQ

if [ "$MODE" = "time" ]; then
  a=$(basename "$(ls "$OUT"/aln/*.fna | head -1)" .fna); n=$(ls "$OUT"/aln/*.fna | wc -l)
  echo "timing anchor: $a   (set size $n, J=$J)"; echo
  for s in cds p12 aa; do
    [ "$s" = aa ] && [ ! -f "$OUT/aa/$a.faa" ] && { echo "  aa   : no peptide alignment, skipped"; continue; }
    t0=$(date +%s.%N); run_one "$s" "$a" >/dev/null 2>&1 || { echo "  $s   : FAILED"; continue; }
    t1=$(date +%s.%N)
    printf "  %-5s: %6.1f core-s/tree  ->  %5.1f core-h  ->  ~%4.0f min wall at -j %s\n" \
      "$s" "$(echo "$t1-$t0" | bc)" "$(echo "($t1-$t0)*$n/3600" | bc -l)" \
      "$(echo "($t1-$t0)*$n/60/$J" | bc -l)" "$J"
  done
  echo; echo "  Above ~30 core-s/tree means something is wrong with the"
  echo "  alignments, not the model. Then: ./02_build.sh cds"
  exit 0
fi

case "$MODE" in
  cds) IN="$OUT/aln"; EXT=fna; TRE="$OUT/tre_cds" ;;
  p12) IN="$OUT/aln"; EXT=fna; TRE="$OUT/tre_p12" ;;
  aa)  IN="$OUT/aa";  EXT=faa; TRE="$OUT/tre_aa"  ;;
  *)   echo "usage: $0 {time|cds|aa|p12}"; exit 1 ;;
esac
n=$(ls "$IN"/*."$EXT" 2>/dev/null | wc -l); [ "$n" -eq 0 ] && { echo "no input in $IN"; exit 1; }
echo "building $MODE : $n anchors at -j $J"; start=$(date +%s)
ls "$IN"/*."$EXT" | xargs -n1 basename | sed "s/\.$EXT\$//" \
  | parallel -j "$J" --bar --joblog "$OUT/log/$MODE.joblog" \
      run_one "$MODE" {} '>' "$OUT/log/$MODE.{}.log" '2>&1'
end=$(date +%s)

echo; echo "======================================================================"
echo "QC — stage 02 ($MODE)"
echo "======================================================================"
echo "  wall time          : $(( (end-start)/60 )) min $(( (end-start)%60 )) s"
echo "  inputs             : $n"
echo "  .treefile produced : $(ls "$TRE"/*.treefile 2>/dev/null | wc -l)"
nz=$(awk 'NR>1 && $7!=0' "$OUT/log/$MODE.joblog" 2>/dev/null | wc -l)
echo "  nonzero exit codes : $nz"
[ "$nz" -gt 0 ] && awk 'NR>1 && $7!=0 {print "    " $NF}' "$OUT/log/$MODE.joblog" | head -5
echo "  empty treefiles    : $(find "$TRE" -name '*.treefile' -empty 2>/dev/null | wc -l)"
echo "  disk used          : $(du -sh "$TRE" 2>/dev/null | cut -f1)"
echo; echo "  example tree:"; head -c 300 "$(ls "$TRE"/*.treefile | head -1)"; echo " ..."
echo; echo "  -> Rscript 03_treeqc.R $MODE"
