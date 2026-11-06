#!/usr/bin/env bash
# launch_expansion_orthofinder.sh -- stage B of the 2026-11 expansion of the holocentricity screen: once stage A1 on
# gpu01 has written all eight final proteomes (logs/expansion/A1.ok, see launch_expansion_gpu01.sh), move the
# 11-species OrthoFinder run aside (once) and run OrthoFinder on all 19 species of config/annotation.csv through Slurm
# (two whole-node jobs, comparative.smk: orthofinder_trees up to the gene trees, then orthofinder from them;
# profiles/slurm/config.yaml, set-resources).
#
# Runs on hpc09 (needs sbatch), started in the background by the setup/fix blocks. Only orthofinder_prep,
# orthofinder_trees and orthofinder may run (--allowed-rules): every species' final proteome is taken as it is on disk,
# so nothing upstream of it is redone for the eleven species annotated earlier, whatever their timestamps. A dry run
# checks that before the real run. Markers in logs/expansion/: B.blocked (dry run planned something else, nothing
# started), B.ok | B.failed.
#
# The 11-species run is kept as results/comparative/orthofinder/out_11sp_<date>/ (with its input and log).
# To start it by hand: nohup setsid bash scripts/launchers/launch_expansion_orthofinder.sh > logs/expansion/B.watch.log 2>&1 < /dev/null &

REPO=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
RUN=$REPO/logs/expansion
SMK=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk
OF=results/comparative/orthofinder
T=$OF/out/Results_drosera/Orthogroups/Orthogroups.tsv
RULES="orthofinder_prep orthofinder_trees orthofinder"

cd "$REPO" || exit 1
export PATH="$SMK/bin:$PATH"   # conda/mamba for the orthofinder env come from the shell that started this
rm -f "$RUN/B.blocked" "$RUN/B.ok" "$RUN/B.failed"
log() { echo "[$(date '+%F %T')] $*"; }

# ---- wait for stage A1 (checked every 10 minutes, for at most 7 days) ----
log "waiting for $RUN/A1.ok"
waited=0
until [ -e "$RUN/A1.ok" ] || [ -e "$RUN/A1.failed" ]; do
  if [ "$waited" -ge 10080 ]; then echo "gave up after 7 days without A1.ok or A1.failed" > "$RUN/B.failed"; log "gave up"; exit 1; fi
  sleep 600; waited=$((waited + 10))
done
if [ -e "$RUN/A1.failed" ]; then
  echo "OrthoFinder not started: stage A1 failed ($(cat "$RUN/A1.failed"))" > "$RUN/B.failed"; log "A1 failed; stopping"; exit 1
fi

# ---- dry run: orthofinder, and before it orthofinder_trees / orthofinder_prep if their outputs are missing or older
#      than their inputs; nothing else ----
snakemake "$T" --profile profiles/slurm --allowed-rules $RULES -n > "$RUN/B.dryrun.txt" 2>&1
rc=$?
planned=$(awk '/^(local)?rule [A-Za-z0-9_]+:/ { r = $2; sub(/:$/, "", r); n[r]++ } END { for (r in n) print r }' "$RUN/B.dryrun.txt" | sort | tr '\n' ' ')
unexpected=$(for r in $planned; do case " $RULES " in *" $r "*) ;; *) echo "$r" ;; esac; done)
if [ "$rc" -ne 0 ] || [ -n "$unexpected" ] || [[ " $planned " != *" orthofinder "* ]]; then
  { echo "OrthoFinder not started: the dry run (exit $rc) planned: ${planned:-nothing}"
    echo "expected: orthofinder, with or without orthofinder_trees / orthofinder_prep. Last lines of $RUN/B.dryrun.txt:"
    tail -n 25 "$RUN/B.dryrun.txt"; } > "$RUN/B.blocked"
  log "blocked: see $RUN/B.blocked"; exit 1
fi

# ---- keep the 11-species run (once: a later restart must not file a failed 19-species run under that name) ----
SNAP=$OF/out_11sp_$(date +%F)
if [ -d "$OF/out" ] && ! ls -d "$OF"/out_11sp_* > /dev/null 2>&1; then
  mv "$OF/out" "$SNAP"
  [ -d "$OF/input" ] && cp -a "$OF/input" "$SNAP/input_11sp"
  [ -s logs/orthofinder/run.log ] && cp -a logs/orthofinder/run.log "$SNAP/run_11sp.log"
  log "11-species run kept as $SNAP"
fi
log "OrthoFinder on $(awk -F, 'NR > 1 && $1 != ""' config/annotation.csv | wc -l) species: submitting through Slurm (${planned})"
# one retry, not the profile's three: a deterministic failure would otherwise cost up to 4 x 4 days
snakemake "$T" --profile profiles/slurm --restart-times 1 --allowed-rules $RULES > "$RUN/B.snakemake.log" 2>&1
rc=$?
n_sp=$(awk -F, 'NR > 1 && $1 != ""' config/annotation.csv | wc -l)
n_col=$([ -s "$T" ] && head -n 1 "$T" | awk -F'\t' '{ print NF - 1 }')
if [ "$rc" -eq 0 ] && [ -s "$T" ] && [ "$n_col" = "$n_sp" ]; then
  echo "OrthoFinder done $(date '+%F %T'): $(($(wc -l < "$T") - 1)) orthogroups, $n_col species in $T" > "$RUN/B.ok"
elif [ "$rc" -eq 0 ] && [ -s "$T" ]; then
  echo "OrthoFinder finished but $T has ${n_col:-?} species columns, not $n_sp: see logs/orthofinder/run.log" > "$RUN/B.failed"
else
  echo "snakemake exit $rc (see $RUN/B.snakemake.log, logs/orthofinder/trees.log and logs/orthofinder/run.log)" > "$RUN/B.failed"
fi
log "stage B finished (snakemake exit $rc)"
