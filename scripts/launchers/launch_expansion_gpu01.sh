#!/usr/bin/env bash
# launch_expansion_gpu01.sh -- stage A of the 2026-11 expansion of the holocentricity screen: annotate the eight
# external genomes (Luzula sylvatica, Juncus effusus, Cuscuta epithymum, Cuscuta campestris, Chionographis japonica,
# Chamaelirium luteum, Morus notabilis, Morus alba) with a LOCAL snakemake on gpu01, which is not in Slurm (see annotation.smk header).
#
#   A1  ANNEVO (GPU) -> proteins -> OMArk vs luca + viridiplantae (species' own taxid) -> keep/drop table
#       -> final GFF3 + proteins. Writes logs/expansion/A1.ok once all eight final proteomes exist, which
#       launch_expansion_orthofinder.sh (hpc09) waits for before starting OrthoFinder.
#   A2  FANTASIA-Lite GO annotation of each final proteome (GPU), while OrthoFinder runs in Slurm.
#
# The genomes were imported (frozen *_chr.fa) on hpc09 by the setup block. Only the rules of this chain may run
# (--allowed-rules), so nothing of the eleven species annotated earlier is touched. GPU jobs share the GPUs through
# --resources gpu=<number of GPUs> plus a per-GPU lock (PICK_GPU, annotation.smk), so two jobs never land on one GPU.
# Markers in logs/expansion/: A.preflight, A1.ok | A1.failed, A2.ok | A2.failed; Snakemake output in A1/A2.snakemake.log.
#
# Started over ssh by the setup block. To start it by hand (e.g. after fixing a failure):
#   ssh gpu01
#   cd /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
#   nohup setsid bash scripts/launchers/launch_expansion_gpu01.sh > logs/expansion/A.log 2>&1 < /dev/null &
# Finished jobs are not redone (--rerun-triggers mtime), so a restart continues where it stopped.

REPO=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
RUN=$REPO/logs/expansion
SMK=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk
NEW="Luzula_sylvatica Juncus_effusus Cuscuta_epithymum Cuscuta_campestris Chionographis_japonica Chamaelirium_luteum Morus_notabilis Morus_alba"
A_RULES="import_external_genome annotate_annevo extract_proteins qc_omark omark_merge omark_classify_genes filter_gff3_by_class final_annotation"

cd "$REPO" || exit 1
mkdir -p "$RUN"
rm -f "$RUN/A.preflight" "$RUN"/A1.ok "$RUN"/A1.failed "$RUN"/A2.ok "$RUN"/A2.failed
# micromamba location and root prefix as found on hpc09 by the setup block (non-login ssh shells may lack them)
[ -s "$RUN/env.sh" ] && . "$RUN/env.sh"
export PATH="$SMK/bin:$PATH"
log() { echo "[$(date '+%F %T')] $*"; }

# ---- preflight: tools, GPUs, cores, memory (the setup block on hpc09 prints A.preflight) ----
missing=""
for c in snakemake python3 micromamba nvidia-smi flock; do command -v "$c" > /dev/null 2>&1 || missing="$missing $c"; done
NGPU=$(nvidia-smi -L 2>/dev/null | grep -c '^GPU')
NPROC=$(nproc)
CORES=$(( NPROC > 68 ? 64 : NPROC - 4 )); [ "$CORES" -ge 8 ] || CORES=$NPROC
MEM=$(awk '/^MemAvailable/ { print int($2 / 1024 * 0.8) }' /proc/meminfo)
if [ -n "$missing" ] || [ "${NGPU:-0}" -lt 1 ] || [ "${MEM:-0}" -lt 70000 ]; then
  echo "FAIL on $(hostname): missing tools:${missing:- none}; GPUs: ${NGPU:-0}; usable memory: ${MEM:-0} MB (ANNEVO and FANTASIA ask 64 GB)" > "$RUN/A.preflight"
  cat "$RUN/A.preflight"; exit 1
fi
{ echo "OK on $(hostname): $NGPU GPU(s), --cores $CORES of $NPROC, --resources mem_mb=$MEM (80% of available)"
  nvidia-smi --query-gpu=index,name,memory.used,memory.total --format=csv,noheader | sed 's/^/   GPU /'; } > "$RUN/A.preflight"
cat "$RUN/A.preflight"

fin() { for s in "$@"; do echo "results/$s/annotation/final/$s.final.proteins.fa"; done; }
SMK_ARGS="--cores $CORES --resources gpu=$NGPU mem_mb=$MEM --rerun-triggers mtime --keep-going --rerun-incomplete --latency-wait 60 --restart-times 1 --printshellcmds"

# ---- A1: annotation; the three largest genomes first (their ANNEVO runs set the end time) ----
BIG=$(for s in $NEW; do f=results/$s/assembly_final/external_collapsed/${s}_chr.fa; [ -s "$f" ] && echo "$(stat -Lc %s "$f") $s"; done \
      | sort -k1,1nr | head -n 3 | awk '{ print $2 }')
PRIO=""; [ -n "$BIG" ] && PRIO="--prioritize $(fin $BIG | tr '\n' ' ')"
log "A1: $NEW (first: $(echo $BIG))"
snakemake $(fin $NEW) $SMK_ARGS $PRIO --allowed-rules $A_RULES > "$RUN/A1.snakemake.log" 2>&1
rc=$?
done_sp=""; miss=""
for s in $NEW; do if [ -s "$(fin $s)" ]; then done_sp="$done_sp $s"; else miss="$miss $s"; fi; done
if [ -z "$miss" ]; then
  echo "snakemake exit $rc; all eight final proteomes written $(date '+%F %T')" > "$RUN/A1.ok"
  log "A1 done (snakemake exit $rc)"
else
  echo "snakemake exit $rc; no final proteome for:$miss  (errors: grep -n Error $RUN/A1.snakemake.log)" > "$RUN/A1.failed"
  log "A1 FAILED for:$miss"
fi

# ---- A2: FANTASIA-Lite for every species that has a final proteome ----
if [ -n "$done_sp" ]; then
  log "A2: FANTASIA for$done_sp"
  snakemake $(for s in $done_sp; do echo "results/$s/annotation/function/fantasia/$s.results.csv"; done) \
    $SMK_ARGS --allowed-rules fantasia > "$RUN/A2.snakemake.log" 2>&1
  rc=$?
  miss2=""; for s in $done_sp; do [ -s "results/$s/annotation/function/fantasia/$s.results.csv" ] || miss2="$miss2 $s"; done
  if [ -z "$miss2" ] && [ -z "$miss" ]; then echo "snakemake exit $rc; FANTASIA done for all eight $(date '+%F %T')" > "$RUN/A2.ok"
  else echo "snakemake exit $rc; no FANTASIA result for:$miss2$miss" > "$RUN/A2.failed"; fi
  log "A2 finished (snakemake exit $rc)"
fi
log "stage A finished"
