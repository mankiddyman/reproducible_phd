#!/usr/bin/env bash
# Stage 01 — deterministic tip rename ('@'->'_' so the treefile matches what we
# wrote) + per-anchor partition files. CDS and P12 share one alignment; P12 is
# just a partition file that omits position 3, so it is free.
set -euo pipefail
BASE="${SUBG_BASE:-$PWD}"
CODON_DIR="$BASE/wgd7/codon_gap"
AA_DIR="${AA_DIR:-}"
OUT="$BASE/trees"
mkdir -p "$OUT"/{aln,aa,part,tre_cds,tre_aa,tre_p12,log}
echo "BASE $BASE"; echo "aa   ${AA_DIR:-<none — AA skipped>}"; echo

n_in=0; n_out=0; n_skip=0
for f in "$CODON_DIR"/*.fna; do
  anchor=$(basename "$f" .fna); n_in=$((n_in+1))
  awk '/^>/ {h=$1; gsub("@","_",h); print h; next} {print}' "$f" > "$OUT/aln/$anchor.fna"
  L=$(awk 'NR>1 && /^>/ {exit} NR>1 {n+=length($0)} END {print n}' "$OUT/aln/$anchor.fna")
  if [ $((L % 3)) -ne 0 ] || [ "$L" -lt 300 ]; then
    n_skip=$((n_skip+1)); rm -f "$OUT/aln/$anchor.fna"; continue
  fi
  printf 'DNA, p1 = 1-%s\\3\nDNA, p2 = 2-%s\\3\nDNA, p3 = 3-%s\\3\n' "$L" "$L" "$L" > "$OUT/part/$anchor.cds.txt"
  printf 'DNA, p1 = 1-%s\\3\nDNA, p2 = 2-%s\\3\n' "$L" "$L" > "$OUT/part/$anchor.p12.txt"
  n_out=$((n_out+1))
done

n_aa=0
if [ -n "$AA_DIR" ] && [ -d "$AA_DIR" ]; then
  for f in "$AA_DIR"/*; do
    case "$f" in *.faa|*.fa|*.fasta|*.aln|*.pep) ;; *) continue ;; esac
    anchor=$(basename "$f"); anchor="${anchor%.*}"
    [ -f "$OUT/aln/$anchor.fna" ] || continue
    awk '/^>/ {h=$1; gsub("@","_",h); print h; next} {print}' "$f" > "$OUT/aa/$anchor.faa"
    n_aa=$((n_aa+1))
  done
fi

echo "======================================================================"
echo "QC — stage 01"
echo "======================================================================"
echo "  codon alignments in      : $n_in"
echo "  prepped (nt)             : $n_out"
echo "  skipped (short / not %3) : $n_skip"
echo "  partition .cds / .p12    : $(ls "$OUT"/part/*.cds.txt 2>/dev/null | wc -l) / $(ls "$OUT"/part/*.p12.txt 2>/dev/null | wc -l)"
echo "  peptide alignments       : $n_aa"
echo; echo "  tip-count check (rename must be 1:1):"
bad=0
for f in "$CODON_DIR"/*.fna; do
  anchor=$(basename "$f" .fna); [ -f "$OUT/aln/$anchor.fna" ] || continue
  a=$(grep -c '^>' "$f"); b=$(grep -c '^>' "$OUT/aln/$anchor.fna")
  c=$(grep '^>' "$OUT/aln/$anchor.fna" | sort -u | wc -l)
  if [ "$a" != "$b" ] || [ "$b" != "$c" ]; then
    echo "    !! $anchor  in=$a out=$b unique=$c"; bad=$((bad+1)); fi
done
[ "$bad" -eq 0 ] && echo "    OK  no tips lost or collided in $n_out alignments" \
                 || { echo "    X   $bad alignments broken by rename"; exit 1; }
echo; echo "  example renamed tips:"
grep '^>' "$(ls "$OUT"/aln/*.fna | head -1)" | head -3 | sed 's/^/    /'
echo; echo "  example partition file:"
sed 's/^/    /' "$(ls "$OUT"/part/*.cds.txt | head -1)"
[ "$n_aa" -eq 0 ] && echo "
  NOTE: no peptide alignments. Set AA_DIR or rebuild with mafft (~15-30 min).
        CDS set is unaffected."
echo "  -> ./02_build.sh time"
