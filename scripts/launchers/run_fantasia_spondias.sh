#!/usr/bin/env bash
#
# run_fantasia_spondias.sh — FANTASIA-Lite functional annotation of the
# collaborator-supplied Spondias tuberosa hap1 + hap2 proteomes.
#
# WHY A SCRIPT AND NOT A SNAKEMAKE RULE
#   Spondias arrives with BOTH assembly and structural annotation already done
#   by collaborators, so it skips the entire annotation track (hifiasm -> decon
#   -> scaffold -> helixer/braker -> omark -> final_annotation) and needs only
#   the terminal FANTASIA step. Wiring that into the DAG needs a new import rule
#   + external_annotations.csv + ruleorder/wildcard_constraints to stop
#   `fantasia` resolving its input through `final_annotation` (which would
#   KeyError on annotation_df for this species). Deferred; this pins provenance.
#
# PROVENANCE
#   Collaborator inputs are NEVER modified in place. Two quirks are normalized
#   by workflow/scripts/normalize_helixer_gtf.py (auto-detecting, see below):
#     1. hap1 GTF seqnames are 'scaffold_N' but its FASTA headers are 'N'.
#        hap2 is already consistent -> the script no-ops on hap2.
#     2. Both GTFs put a functional DESCRIPTION in transcript_id, which is
#        non-unique (hap1: 24151 genes -> 15971 distinct transcript_id strings;
#        hap2: 23811 -> 14053). True transcript IDs are recovered from the ID
#        attribute (mRNA lines) / Parent attribute (exon, CDS, UTR lines).
#   gene_id is NEVER touched, so FANTASIA output joins back to collaborator gene
#   names (hap1_scaffold_1_003111 / hap2_1_003098) and to the CO_smk landscape.
#
# USAGE  (must be on a GPU node — ProtT5 runs on cuda)
#   ssh gpu01
#   cd /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
#   bash scripts/launchers/run_fantasia_spondias.sh          # both haps
#   bash scripts/launchers/run_fantasia_spondias.sh hap1     # one hap

set -euo pipefail

GDD=/netscratch/dep_mercier/grp_marques/marques/Spondias_tuberosa/google_drive_download
AGAT_ENV=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/agat_env
FANTASIA_ENV=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/fantasia-lite
FANTASIA_PIPE=/netscratch/dep_mercier/grp_marques/Aaryan/methods/FANTASIA-Lite/src/fantasia_pipeline.py
TMPBASE=/netscratch/dep_mercier/grp_marques/Aaryan/tmp/fantasia
REPO_ROOT=$(pwd)

HAPS=("${@:-hap1 hap2}")
[ $# -gt 0 ] && HAPS=("$@") || HAPS=(hap1 hap2)

# --- preflight -------------------------------------------------------------
[ -f "$FANTASIA_PIPE" ] || { echo "FATAL: FANTASIA pipeline not found: $FANTASIA_PIPE"; exit 1; }
nvidia-smi -L >/dev/null 2>&1 || { echo "FATAL: no GPU visible. FANTASIA needs cuda -> ssh gpu01 first."; exit 1; }

GPU_ID=$(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits \
         | sort -t, -k2 -n | head -1 | cut -d, -f1 | tr -d ' ')
export CUDA_VISIBLE_DEVICES="$GPU_ID"
echo "Using GPU $GPU_ID"

for H in "${HAPS[@]}"; do
    SP="Spondias_tuberosa_${H}"
    FA="$GDD/${H}scaffold.FINAL.fasta"
    GTF_SRC="$GDD/${H}scaffold.FINAL_helixer.gtf"

    ANNOT_DIR="$REPO_ROOT/results/$SP/annotation"
    GTF_NORM="$ANNOT_DIR/external/${SP}.normalized.gtf"
    FAA="$ANNOT_DIR/final/${SP}.final.proteins.fa"
    OUT="$ANNOT_DIR/function/fantasia"
    LOG="$REPO_ROOT/logs/fantasia/${SP}.log"

    mkdir -p "$(dirname "$GTF_NORM")" "$(dirname "$FAA")" "$OUT/topgo" "$(dirname "$LOG")"

    echo "===== $SP ====="
    [ -f "$FA" ]      || { echo "FATAL: missing $FA"; exit 1; }
    [ -f "$GTF_SRC" ] || { echo "FATAL: missing $GTF_SRC"; exit 1; }
    [ -f "$FA.fai" ]  || samtools faidx "$FA"

    # 1. normalize collaborator GTF (idempotent; skip if present)
    if [ ! -s "$GTF_NORM" ]; then
        echo "--- normalizing GTF ---"
        python3 "$REPO_ROOT/workflow/scripts/normalize_helixer_gtf.py" \
            --gtf "$GTF_SRC" --fasta "$FA" --out "$GTF_NORM"
    else
        echo "--- normalized GTF exists, skipping ---"
    fi

    # 2. extract proteome (idempotent; skip if present)
    if [ ! -s "$FAA" ]; then
        echo "--- extracting proteome (AGAT) ---"
        micromamba run -p "$AGAT_ENV" agat_sp_extract_sequences.pl \
            -g "$GTF_NORM" -f "$FA" -p -o "$FAA" > "$REPO_ROOT/logs/fantasia/${SP}.agat.log" 2>&1
    else
        echo "--- proteome exists, skipping ---"
    fi
    echo "proteins: $(grep -c '^>' "$FAA")"

    # 3. FANTASIA (chunk I/O on netscratch, NOT local /scratch)
    mkdir -p "$TMPBASE"
    SCRATCH=$(mktemp -d "$TMPBASE/${SP}_XXXXXX")
    trap 'rm -rf "$SCRATCH"' EXIT

    echo "--- FANTASIA (GPU $CUDA_VISIBLE_DEVICES) -> $LOG ---"
    echo "=== FANTASIA $SP (GPU $CUDA_VISIBLE_DEVICES) $(date) ===" >> "$LOG"
    micromamba run -p "$FANTASIA_ENV" python "$FANTASIA_PIPE" \
        "$FAA" \
        --device cuda \
        --results-csv     "$OUT/${SP}.results.csv" \
        --topgo-dir       "$OUT/topgo" \
        --embeddings-npz  "$OUT/${SP}.query_embeddings.npz" \
        --config-yaml     "$OUT/${SP}.fantasia_config.yaml" \
        --failure-report  "$OUT/${SP}.failed_sequences.csv" \
        --chunk-dir          "$SCRATCH/fasta_chunks" \
        --chunk-embed-dir    "$SCRATCH/chunk_embeddings" \
        --chunk-results-dir  "$SCRATCH/chunk_results" \
        --chunk-config-dir   "$SCRATCH/chunk_configs" \
        --chunk-failure-dir  "$SCRATCH/chunk_failures" \
        >> "$LOG" 2>&1

    rm -rf "$SCRATCH"; trap - EXIT
    echo "done: $OUT/${SP}.results.csv"
    wc -l "$OUT/${SP}.results.csv" || true
done

echo "=== all done ==="
