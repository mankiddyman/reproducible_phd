#!/bin/bash
set -euo pipefail
ROOT=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
OUT=$ROOT/results/comparative/subgenome_prototype/cds
GSD=$ROOT/results/comparative/genespace/results
mkdir -p "$OUT"

tail -n +2 "$ROOT/config/genespace.csv" | tr -d '\r' | while IFS=, read -r sp stage regex gsname ploidy; do
  [ -z "${sp:-}" ] && continue
  GEN="$ROOT/results/$sp/assembly_final/$stage/${sp}_chr.fa"
  SUB="$ROOT/results/comparative/genespace/gff/${sp}.gff3"      # the SUBSET gff GENESPACE used
  FULL="$ROOT/results/$sp/annotation/final/${sp}.final.gff3"
  PEP="$ROOT/results/comparative/genespace/peptide/${sp}.fa"

  # ASSERTION 1 - genome exists at the stage named in config
  [ -s "$GEN" ] || { echo "MISSING genome for $sp: $GEN" >&2; exit 1; }
  [ -s "$PEP" ] || { echo "MISSING genespace peptide for $sp" >&2; exit 1; }

  # prefer the stored subset gff; else rebuild it with the SAME regex
  if [ -s "$SUB" ]; then GFF="$SUB"
  else
    [ -s "$FULL" ] || { echo "MISSING gff for $sp" >&2; exit 1; }
    GFF="$OUT/${sp}.subset.gff3"
    awk -F'\t' -v re="$regex" '/^#/{print; next} $1 ~ re' "$FULL" > "$GFF"
  fi

  # ASSERTION 2 - genome contains every chromosome combBed has for this species
  [ -s "$GEN.fai" ] || samtools faidx "$GEN"
  comm -13 <(cut -f1 "$GEN.fai" | sort -u) \
           <(awk -F'\t' -v s="$sp" 'NR>1 && $6==s {print $1}' "$GSD/combBed.txt" | sort -u) > "$OUT/.miss"
  if [ -s "$OUT/.miss" ]; then
    echo "GENOME MISMATCH for $sp ($stage): combBed chromosomes absent from $GEN" >&2
    head -5 "$OUT/.miss" >&2; exit 1
  fi

  gffread "$GFF" -g "$GEN" -x "$OUT/${sp}.cds.fa" 2>/dev/null

  # ASSERTION 3 - one CDS per GENESPACE peptide, same IDs
  NP=$(grep -c '^>' "$PEP"); NC=$(grep -c '^>' "$OUT/${sp}.cds.fa")
  SHARED=$(comm -12 <(grep '^>' "$PEP" | sed 's/^>//;s/[[:space:]].*//' | sort) \
                    <(grep '^>' "$OUT/${sp}.cds.fa" | sed 's/^>//;s/[[:space:]].*//' | sort) | wc -l)
  printf '%-22s stage=%-19s pep=%-6s cds=%-6s shared=%-6s\n' "$sp" "$stage" "$NP" "$NC" "$SHARED"
  [ "$SHARED" -ge $((NP * 95 / 100)) ] || { echo "  ID MISMATCH for $sp" >&2; exit 1; }

  seqkit replace -p '^' -r "${sp}@" "$OUT/${sp}.cds.fa" > "$OUT/${sp}.cds.tagged.fa"
done
rm -f "$OUT/.miss"
cat "$OUT"/*.cds.tagged.fa > "$OUT/all_cds_tagged.fa"
echo; seqkit stats "$OUT/all_cds_tagged.fa"
