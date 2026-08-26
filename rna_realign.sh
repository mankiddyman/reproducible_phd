#!/usr/bin/env bash
# Realign RNA-seq per species -> per-exon coverage vs final.gff3.
# Standalone (NOT a snakemake rule): ad-hoc annotation-trust assessment.
# usage: bash rna_realign.sh <species>
set -euo pipefail

SP="$1"
RP=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
ENV=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk
OUT="$RP/results/$SP/rna_support"
THREADS=24

cd "$RP"
export PATH="$ENV/bin:$PATH"
export MAMBA_ROOT_PREFIX=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba
mkdir -p "$OUT"
LOG="$OUT/run.log"; : > "$LOG"
echo "=== $SP  $(date) ===" | tee -a "$LOG"

# --- inputs ---
GFF="$RP/results/$SP/annotation/final/$SP.final.gff3"
STAGE=$(awk -F, -v s="$SP" '$1==s{print $7}' config/annotation.csv)
GEN="$RP/results/$SP/assembly_final/$STAGE/${SP}_chr.fa"
RNADIR=$(awk -F, -v s="$SP" '$1==s{print $3}' config/annotation.csv)
for f in "$GFF" "$GEN"; do [ -f "$f" ] || { echo "MISSING $f" | tee -a "$LOG"; exit 1; }; done
[ -d "$RNADIR" ] || { echo "MISSING rna dir $RNADIR" | tee -a "$LOG"; exit 1; }

R1=$(ls "$RNADIR"/*_R1_*.fastq.gz | sort | paste -sd, -)
R2=$(echo "$R1" | sed 's/_R1_/_R2_/g')
echo "genome:  $GEN"      | tee -a "$LOG"
echo "R1:      $R1"       | tee -a "$LOG"

command -v hisat2 >/dev/null || { echo "hisat2 not on PATH" | tee -a "$LOG"; exit 1; }
command -v samtools >/dev/null || { echo "samtools not on PATH" | tee -a "$LOG"; exit 1; }

# --- 1. index ---
IDX="$OUT/idx"
if [ ! -f "$IDX.1.ht2" ] && [ ! -f "$IDX.1.ht2l" ]; then
  echo "[$(date +%T)] hisat2-build" | tee -a "$LOG"
  hisat2-build -p $THREADS "$GEN" "$IDX" >> "$LOG" 2>&1
fi

# --- 2. align (CSI index: chromosomes >512 Mb) ---
BAM="$OUT/$SP.rna.sorted.bam"
if [ ! -f "$BAM" ]; then
  echo "[$(date +%T)] hisat2 align" | tee -a "$LOG"
  hisat2 -p $THREADS -x "$IDX" -1 "$R1" -2 "$R2" \
      --dta --new-summary --summary-file "$OUT/hisat2.summary" 2>> "$LOG" \
    | samtools sort -@ 8 -m 4G -o "$BAM" - 2>> "$LOG"
  samtools index -c -@ 8 "$BAM" 2>> "$LOG"
fi

# --- 3. per-exon coverage ---
echo "[$(date +%T)] per-exon coverage" | tee -a "$LOG"
awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && $3=="exon" {
    id=""; if (match($9,/Parent=[^;]+/)) id=substr($9,RSTART+7,RLENGTH-7)
    print $1, $4-1, $5, id, ".", $7 }' "$GFF" | sort -k1,1 -k2,2n > "$OUT/exons.bed"
echo "exons: $(wc -l < "$OUT/exons.bed")" | tee -a "$LOG"

samtools bedcov "$OUT/exons.bed" "$BAM" > "$OUT/exon_bedcov.tsv" 2>> "$LOG"
samtools view -b -q 10 "$BAM" 2>/dev/null | bedtools coverage -a "$OUT/exons.bed" -b stdin \
    > "$OUT/exon_breadth.tsv" 2>> "$LOG" || echo "(bedtools missing - breadth skipped)" | tee -a "$LOG"

echo "[$(date +%T)] DONE $SP" | tee -a "$LOG"
