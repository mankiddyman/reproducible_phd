#!/bin/bash
set -euo pipefail
cd /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome
GENOME=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/Dionaea_muscipula/assembly_final/external_collapsed/Dionaea_muscipula_chr.fa
PEP=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/genespace/peptide/Dionaea_muscipula.fa
THREADS=8
PAR=12

while read -r c; do
  [ -s blast/db/$c.nsq ] && continue
  samtools faidx "$GENOME" "$c" > blast/chr/$c.fa
  makeblastdb -in blast/chr/$c.fa -dbtype nucl -out blast/db/$c -title "$c" >/dev/null
done < blast/chroms.txt

for f in blast/prot/*.ids; do
  t=$(basename "$f" .ids)
  [ -s blast/prot/$t.fa ] || seqkit grep -f "$f" "$PEP" > blast/prot/$t.fa
done

tail -n +2 blast/jobs.csv | awk -F, '{print $1","$4}' | \
xargs -P $PAR -I{} bash -c '
  IFS=, read -r tag tgt <<< "{}"
  [ -s blast/out/$tag.tsv ] && exit 0
  tblastn -query blast/prot/$tag.fa -db blast/db/$tgt \
    -outfmt "6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore" \
    -evalue 1e-5 -max_target_seqs 50 -num_threads '"$THREADS"' \
    -out blast/out/$tag.tsv 2> blast/out/$tag.err
  echo "done $tag -> $tgt"
'
