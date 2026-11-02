#!/usr/bin/env bash
# tpxl6_te_loci.sh -- repeat annotation (CARP) in and around every TPXL6 copy
#
# For every species with a TPXL6 copy (tpx2/miniprot/tpxl6_miniprot_loci.tsv, assigned TPXL6) and a finished CARP
# run (results/{species}/repeats/carp/Repeat_Annotation_Unified.gff3), keeps the CARP Level-1 features (no Parent:
# one classification per base) that overlap a TPXL6 copy or lie within FLANK bp of it.
# A species is skipped, and says why, when CARP has not finished or ran on a different assembly from the one
# miniprot searched (genome_fasta in results/{species}/repeats/carp_run/carp_config.yaml).
# Output: results/comparative/holocentricity/tpx2/te_loci/{species}.te.gff and status.tsv
# usage: bash tpxl6_te_loci.sh [REPO] [FLANK_BP]
REPO=${1:-/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd}
FL=${2:-5000}
MP=$REPO/results/comparative/holocentricity/tpx2/miniprot
OUT=$REPO/results/comparative/holocentricity/tpx2/te_loci
Q6=regia_g4881.t1
cd "$REPO" || exit 1
[ -s "$MP/tpxl6_miniprot_loci.tsv" ] || { echo "!! $MP/tpxl6_miniprot_loci.tsv missing: run tpxl6_miniprot.sh first"; exit 1; }
mkdir -p "$OUT"
echo -e "species\tcarp\tnote" > "$OUT/status.tsv"
loci() {  # contig, start, end of every TPXL6 copy of species $1
  awk -F'\t' -v sp="$1" -v q6="$Q6" 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
    $c["species"] == sp && $c["query"] == q6 && $c["assigned"] == "TPXL6" { print $c["contig"] "\t" $c["start"] "\t" $c["end"] }' \
    "$MP/tpxl6_miniprot_loci.tsv"
}
for sp in $(awk -F'\t' -v q6="$Q6" 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
              $c["query"] == q6 && $c["assigned"] == "TPXL6" { print $c["species"] }' "$MP/tpxl6_miniprot_loci.tsv" | sort -u); do
  G=results/$sp/repeats/carp/Repeat_Annotation_Unified.gff3
  if [ ! -s "$G" ]; then
    echo -e "$sp\tno\tno CARP annotation yet" >> "$OUT/status.tsv"; echo "== $sp: no CARP annotation yet"; continue
  fi
  mp_genome=$(readlink -f "$(awk -F'\t' -v sp="$sp" '$1 == sp { print $2 }' "$MP/genomes.tsv")" 2>/dev/null)
  carp_raw=$(awk '/^genome_fasta:/ { print $2 }' "results/$sp/repeats/carp_run/carp_config.yaml" 2>/dev/null)
  carp_genome=${carp_raw:+$(readlink -f "$carp_raw" 2>/dev/null || echo "$carp_raw")}
  if [ -n "$carp_genome" ] && [ "$carp_genome" != "$mp_genome" ]; then
    echo -e "$sp\tno\tCARP ran on $carp_genome, miniprot on $mp_genome" >> "$OUT/status.tsv"
    echo "== $sp: skipped, CARP ran on $carp_genome but miniprot on $mp_genome"; continue
  fi
  loci "$sp" | awk -v fl="$FL" -F'\t' '{ print $1 "\t" $2 - fl "\t" $3 + fl }' > "$OUT/$sp.regions"
  awk -F'\t' 'NR == FNR { n++; C[n] = $1; L[n] = $2; H[n] = $3; want[$1] = 1; next }
              /^#/ || !($1 in want) || $9 ~ /(^|;)Parent=/ { next }
              { for (i = 1; i <= n; i++) if ($1 == C[i] && $4 <= H[i] && $5 >= L[i]) { print; break } }' \
    "$OUT/$sp.regions" "$G" > "$OUT/$sp.te.gff"
  note="$(wc -l < "$OUT/$sp.te.gff") CARP features within $FL bp of $(wc -l < "$OUT/$sp.regions") TPXL6 copies"
  [ -n "$carp_genome" ] || note="$note (assembly not checked: no carp_config.yaml)"
  echo -e "$sp\tyes\t$note" >> "$OUT/status.tsv"; echo "== $sp: $note"
done
