#!/usr/bin/env bash
# =====================================================================
# Gate 1: splitting-artifact verification for the 7 genome-biology
# P/A holocentricity candidate OGs (parked hypothesis, post-Copenhagen).
#
# Logic: a "present in group A, absent in group B" OG may be a real
# gain/loss OR an OrthoFinder orthogroup-splitting artifact (the gene
# exists in group B but landed in a different OG). BLASTp the EXISTING
# proteins against the ABSENT group's proteomes:
#   strong hit in absent group => homolog exists => SPLITTING ARTIFACT
#   no/weak hit                 => absence may be real => keep candidate
#
# Bidirectional:
#   6 holo-gain OGs: query = HOLO proteins -> search MONO proteomes
#   1 holo-loss OG (OG0012850): query = MONO proteins -> search HOLO proteomes
#
# Run from repo root. Needs: seqkit (module), makeblastdb+blastp (/usr/bin).
# =====================================================================
set -euo pipefail

REPO="/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd"
WORK="$REPO/results/comparative/holocentricity/blast_verify"
MEMBERS="$REPO/results/comparative/holocentricity/pa_strict_chrom_members.tsv"
THREADS=8

HOLO=(Drosera_paradoxa Drosera_regia Drosera_roseana)
MONO=(Drosera_aliciae Drosera_binata Drosera_capensis Drosera_filiformis Drosera_tokaiensis)

# strong-hit thresholds (homolog considered present in the absent group)
MIN_PID=50      # % identity
MAX_EVAL=1e-20  # e-value
MIN_QCOV=50     # % query coverage

mkdir -p "$WORK"; cd "$WORK"
proteome() { echo "$REPO/results/$1/annotation/final/$1.final.proteins.fa"; }

echo "### STEP 1: extract query gene ids per OG (split by direction) ###"
# member cols: corner Orthogroup species n_genes genes(comma-sep)
awk -F'\t' 'NR>1{n=split($5,g,","); for(i=1;i<=n;i++) print $2"\t"$1"\t"$3"\t"g[i]}' \
  "$MEMBERS" > og_gene_long.tsv
echo "  query gene rows: $(wc -l < og_gene_long.tsv)"

echo "### STEP 2: build a combined proteome with bare gene-id headers ###"
# headers ">g123.t1 gene=g123 ..."; rewrite to ">SPECIES|g123" so we can both
# extract queries by id AND know which species each query came from.
: > all_proteomes_tagged.faa
for sp in "${HOLO[@]}" "${MONO[@]}"; do
  seqkit replace -p '.*gene=(\S+).*' -r "${sp}|\${1}" "$(proteome "$sp")" \
    >> all_proteomes_tagged.faa
done
echo "  combined tagged proteome seqs: $(grep -c '^>' all_proteomes_tagged.faa)"

echo "### STEP 3: pull QUERY sequences ###"
# holo-gain queries = HOLO species genes in holo-gain OGs
awk -F'\t' '$2=="holo-gain"{print $3"|"$4}' og_gene_long.tsv \
  | grep -E "^(Drosera_paradoxa|Drosera_regia|Drosera_roseana)\|" | sort -u > q_gain.tagids
# holo-loss queries = MONO species genes in the holo-loss OG
awk -F'\t' '$2=="holo-loss"{print $3"|"$4}' og_gene_long.tsv \
  | grep -E "^(Drosera_aliciae|Drosera_binata|Drosera_capensis|Drosera_filiformis|Drosera_tokaiensis)\|" | sort -u > q_loss.tagids
echo "  holo-gain query ids: $(wc -l < q_gain.tagids); holo-loss query ids: $(wc -l < q_loss.tagids)"

seqkit grep -f q_gain.tagids all_proteomes_tagged.faa > query_gain.faa
seqkit grep -f q_loss.tagids all_proteomes_tagged.faa > query_loss.faa
echo "  query_gain.faa seqs: $(grep -c '^>' query_gain.faa) (expect $(wc -l < q_gain.tagids))"
echo "  query_loss.faa seqs: $(grep -c '^>' query_loss.faa) (expect $(wc -l < q_loss.tagids))"

echo "### STEP 4: build BLAST dbs of the ABSENT-group proteomes ###"
mkdir -p db
# holo-gain absent group = MONO ; holo-loss absent group = HOLO
for sp in "${MONO[@]}" "${HOLO[@]}"; do
  seqkit replace -p '.*gene=(\S+).*' -r "${sp}|\${1}" "$(proteome "$sp")" > db/${sp}.faa
  makeblastdb -in db/${sp}.faa -dbtype prot -out db/${sp} >/dev/null 2>&1
done
echo "  built $(ls db/*.pin 2>/dev/null | wc -l) blast dbs"

echo "### STEP 5: blastp query vs each absent-group proteome ###"
# outfmt: query, subject, %id, alnlen, evalue, bitscore, qcovs
BFMT="6 qseqid sseqid pident length evalue bitscore qcovs"
: > blast_gain.tsv
for sp in "${MONO[@]}"; do
  blastp -query query_gain.faa -db db/$sp -outfmt "$BFMT" \
    -evalue 1 -max_target_seqs 5 -num_threads $THREADS \
    | awk -v s="$sp" 'BEGIN{OFS="\t"}{print s,$0}' >> blast_gain.tsv
done
: > blast_loss.tsv
for sp in "${HOLO[@]}"; do
  blastp -query query_loss.faa -db db/$sp -outfmt "$BFMT" \
    -evalue 1 -max_target_seqs 5 -num_threads $THREADS \
    | awk -v s="$sp" 'BEGIN{OFS="\t"}{print s,$0}' >> blast_loss.tsv
done
echo "  blast_gain hits: $(wc -l < blast_gain.tsv); blast_loss hits: $(wc -l < blast_loss.tsv)"

echo "### STEP 6: classify per OG (artifact vs possibly-real) ###"
# map each query gene -> its OG
awk -F'\t' 'NR>1{n=split($5,g,","); for(i=1;i<=n;i++) print $3"|"g[i]"\t"$2}' \
  "$MEMBERS" > tagid2og.tsv

# for each (OG), does ANY query gene get a STRONG hit in the absent group?
# strong = pident>=MIN_PID & evalue<=MAX_EVAL & qcovs>=MIN_QCOV
classify() {
  local blast="$1" corner="$2"
  # cols of blast file: db_species qseqid sseqid pident length evalue bitscore qcovs
  awk -F'\t' -v p=$MIN_PID -v e=$MAX_EVAL -v q=$MIN_QCOV '
    { pid=$4; ev=$6+0; qc=$8;
      if (pid>=p && ev<=e && qc>=q) print $2 }   # $2 = qseqid (SPECIES|gene)
  ' "$blast" | sort -u > /tmp/strong_$$.ids
  # join strong query ids -> OG, count strong hits per OG
  join -t$'\t' -1 1 -2 1 \
    <(sort -k1,1 tagid2og.tsv) \
    <(sort /tmp/strong_$$.ids) 2>/dev/null \
    | awk -F'\t' '{c[$2]++} END{for(o in c) print o"\t"c[o]}' | sort > /tmp/strong_by_og_$$.tsv
  rm -f /tmp/strong_$$.ids
}

echo
echo "================= RESULT ================="
echo "Per-OG verdict (strong hit in ABSENT group => SPLITTING ARTIFACT):"
echo "thresholds: pident>=$MIN_PID, evalue<=$MAX_EVAL, qcov>=$MIN_QCOV"
echo "------------------------------------------"
# all candidate OGs + corner
awk -F'\t' 'NR>1{print $2"\t"$1}' "$MEMBERS" | sort -u > all_ogs.tsv

classify blast_gain.tsv holo-gain; cp /tmp/strong_by_og_$$.tsv strong_gain.tsv 2>/dev/null || : > strong_gain.tsv
classify blast_loss.tsv holo-loss; cp /tmp/strong_by_og_$$.tsv strong_loss.tsv 2>/dev/null || : > strong_loss.tsv
rm -f /tmp/strong_by_og_$$.tsv

cat strong_gain.tsv strong_loss.tsv > strong_all.tsv 2>/dev/null || : > strong_all.tsv
while IFS=$'\t' read -r og corner; do
  nstrong=$(awk -F'\t' -v o="$og" '$1==o{print $2}' strong_all.tsv)
  nstrong=${nstrong:-0}
  if [ "$nstrong" -gt 0 ]; then
    verdict="ARTIFACT (homolog in absent group; $nstrong query genes hit strongly)"
  else
    verdict="possibly-real (no strong hit in absent group)"
  fi
  printf "%-12s %-10s %s\n" "$og" "$corner" "$verdict"
done < all_ogs.tsv

echo "------------------------------------------"
echo "Full hit tables: blast_gain.tsv, blast_loss.tsv"
echo "Inspect any OG's best hits e.g.:"
echo "  grep <gene_id> blast_gain.tsv | sort -t\$'\t' -k5,5 -nr | head"
