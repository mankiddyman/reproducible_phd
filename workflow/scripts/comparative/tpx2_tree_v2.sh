#!/usr/bin/env bash
# tpx2_tree_v2.sh -- the TPX2-family gene tree, rebuilt with every step recorded and with bootstrap support
#
# Why a rebuild: the June tree (tpx2/tpx2_family.nwk) carries 0-1 support values, not bootstrap percentages,
# and no script in git records how it was built.
#
# What it does (all outputs in results/comparative/holocentricity/tpx2/v2/):
#   1. inputs: the same 237 proteins (tpx2/tpx2_family.faa) plus the D. binata TPXL6 reading frame, taken from
#      tpx2/tpx2_short.faa (where its tip is called binata_PSEUDOGENE) and renamed binata_TPXL6_ORF.
#      Headers cut to the ID, terminal stops removed, internal stops (none expected) turned into X and reported.
#   2. MAFFT L-INS-i: suited to one alignable domain (here the TPX2 C-terminal domain) with variable flanks
#   3. trimAl -automated1: drops gap-rich, poorly aligned columns; reports how many columns survive and which
#      proteins keep less than 30% of the trimmed alignment (their placement is the least reliable)
#   4. IQ-TREE: ModelFinder picks the substitution model; 1000 ultrafast bootstraps (refined with --bnni)
#      and 1000 SH-aLRT replicates, so every node carries "SH-aLRT/UFBoot".
#      Usual reading: a split is well supported when SH-aLRT >= 80 and UFBoot >= 95.
#   5. PROVENANCE.txt: tool versions, input checksums, counts, chosen model, commands
#
# usage: bash tpx2_tree_v2.sh [REPO] [THREADS]      (Slurm job: see the commit that added this script)
REPO=${1:-/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd}
THREADS=${2:-32}
ENV=${PHYLO_ENV:-/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/phylo}
B=$ENV/bin
IN=$REPO/results/comparative/holocentricity/tpx2
OUT=$IN/v2
die() { echo "!! $*" >&2; exit 1; }
for t in mafft trimal iqtree3 seqkit; do [ -x "$B/$t" ] || die "$B/$t missing: create the env from envs/phylo.yaml"; done
for f in "$IN/tpx2_family.faa" "$IN/tpx2_short.faa"; do [ -s "$f" ] || die "missing input $f"; done
mkdir -p "$OUT" && cd "$OUT" || die "cannot enter $OUT"
T0=$(date +%s)

# FASTA -> IDs only, one line per sequence, terminal '*' removed, internal '*' -> X (reported on stderr)
clean_fa() {
  awk 'function out() { if (id == "") return; s = seq; sub(/\*+$/, "", s); n = gsub(/\*/, "X", s)
                        if (n) printf "   %s: %d internal stop(s) turned into X\n", id, n > "/dev/stderr"
                        printf ">%s\n%s\n", id, s }
       /^>/ { out(); id = substr($1, 2); seq = ""; next }
       { gsub(/[ \t\r]/, ""); seq = seq $0 }
       END { out() }' "$1"
}

echo "== 1. inputs   $(date '+%F %T')"
clean_fa "$IN/tpx2_family.faa" > family.faa
clean_fa "$IN/tpx2_short.faa" | awk '/^>/ { keep = ($0 ~ /PSEUDOGENE/) } keep' > binata_raw.faa
nb=$(grep -c '^>' binata_raw.faa)
[ "$nb" -eq 1 ] || die "expected exactly 1 PSEUDOGENE record in tpx2_short.faa, found $nb"
awk 'NR == 1 { print ">binata_TPXL6_ORF"; next } { print }' binata_raw.faa > binata_TPXL6_ORF.faa
grep -q '^>binata_TPXL6_ORF$' family.faa && die "binata_TPXL6_ORF already in tpx2_family.faa"
cat family.faa binata_TPXL6_ORF.faa > tpx2_v2.faa
dup=$(grep '^>' tpx2_v2.faa | sort | uniq -d | head -3)
[ -z "$dup" ] || die "duplicate IDs: $dup"
NSEQ=$(grep -c '^>' tpx2_v2.faa)
echo "   $(grep -c '^>' family.faa) family proteins + binata_TPXL6_ORF = $NSEQ sequences"
awk '/^>/ { next } { L = length($0); x = gsub(/X/, "X"); m = substr($0, 1, 1)
       printf "   binata_TPXL6_ORF: %d aa, starts with %s, %d X (internal stops)\n", L, m, x }' binata_TPXL6_ORF.faa

echo "== 2. MAFFT L-INS-i   $(date '+%F %T')"
"$B/mafft" --localpair --maxiterate 1000 --thread "$THREADS" tpx2_v2.faa > tpx2_v2.mafft.faa 2> mafft.log \
  || die "mafft failed (mafft.log)"
[ "$(grep -c '^>' tpx2_v2.mafft.faa)" -eq "$NSEQ" ] || die "alignment does not hold $NSEQ sequences"
ALN=$(awk '/^>/ { if (n++) exit; next } { L += length($0) } END { print L }' tpx2_v2.mafft.faa)
echo "   $ALN alignment columns"

echo "== 3. trimAl -automated1   $(date '+%F %T')"
"$B/trimal" -in tpx2_v2.mafft.faa -out tpx2_v2.trim.faa -fasta -automated1 -colnumbering > trimal_kept_columns.txt 2> trimal.log \
  || die "trimal failed (trimal.log)"
TRIM=$(awk '/^>/ { if (n++) exit; next } { L += length($0) } END { print L }' tpx2_v2.trim.faa)
echo "   $TRIM of $ALN columns kept"
awk -v T="$TRIM" '/^>/ { if (id != "") print id "\t" r; id = substr($1, 2); r = 0; next }
                  { L = length($0); g = gsub(/[-X?]/, ""); r += L - g }
                  END { print id "\t" r }' tpx2_v2.trim.faa \
  | awk -v T="$TRIM" -F'\t' 'BEGIN { OFS = "\t"; print "id", "residues_in_trimmed_alignment", "fraction" }
                            { print $1, $2, sprintf("%.2f", $2 / T) }' > trimmed_occupancy.tsv
LOW=$(awk -F'\t' 'NR > 1 && $3 < 0.30' trimmed_occupancy.tsv | wc -l)
echo "   $LOW proteins keep < 30% of the trimmed columns (listed in trimmed_occupancy.tsv):"
awk -F'\t' 'NR > 1 && $3 < 0.30 { printf "     %s %s\n", $1, $3 }' trimmed_occupancy.tsv | head -20
awk -F'\t' '$1 == "binata_TPXL6_ORF" { printf "   binata_TPXL6_ORF keeps %s of the trimmed columns\n", $3 }' trimmed_occupancy.tsv

echo "== 4. IQ-TREE (ModelFinder, 1000 UFBoot + bnni, 1000 SH-aLRT)   $(date '+%F %T')"
"$B/iqtree3" -s tpx2_v2.trim.faa -m MFP -B 1000 --bnni --alrt 1000 -T AUTO --threads-max "$THREADS" \
  --seed 20261030 --prefix tpx2_v2 --redo > iqtree.stdout 2>&1 || die "iqtree failed (iqtree.stdout, tpx2_v2.log)"
[ -s tpx2_v2.treefile ] || die "no treefile"
MODEL=$(grep -m1 'Best-fit model' tpx2_v2.iqtree | sed 's/.*: //')
echo "   best model: $MODEL"

echo "== 5. provenance   $(date '+%F %T')"
{
  echo "tpx2_tree_v2.sh run $(date '+%F %T') on $(hostname), $THREADS threads, $(( ($(date +%s) - T0) / 60 )) min"
  echo "git commit: $(git -C "$REPO" rev-parse --short HEAD 2>/dev/null)"
  echo "inputs:"; md5sum "$IN/tpx2_family.faa" "$IN/tpx2_short.faa" | sed 's/^/  /'
  echo "sequences: $NSEQ ($(( NSEQ - 1 )) family + binata_TPXL6_ORF from the PSEUDOGENE record of tpx2_short.faa)"
  echo "alignment: MAFFT $("$B/mafft" --version 2>&1 | head -1), --localpair --maxiterate 1000 -> $ALN columns"
  echo "trimming: $("$B/trimal" --version 2>&1 | grep -m1 -i trimal), -automated1 -> $TRIM columns; $LOW proteins < 30% occupancy"
  echo "tree: $("$B/iqtree3" --version 2>&1 | head -1)"
  echo "  -m MFP -B 1000 --bnni --alrt 1000 -T AUTO --threads-max $THREADS --seed 20261030"
  echo "  best-fit model: $MODEL"
  echo "  node labels: SH-aLRT/UFBoot"
  echo "env: $ENV"
  ls "$ENV/conda-meta" 2>/dev/null | grep -E '^(mafft|trimal|iqtree|seqkit)-[0-9]' | sed 's/\.json$//; s/^/  package: /'
} > PROVENANCE.txt
cat PROVENANCE.txt
echo "== done   $(date '+%F %T')   tree: $OUT/tpx2_v2.treefile"
