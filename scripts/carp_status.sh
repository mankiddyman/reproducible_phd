#!/usr/bin/env bash
# carp_status.sh -- where every species' CARP repeat annotation stands, one line each.
#   DONE     published (MANIFEST.txt written by rule publish_repeats)
#   RUNNING  a Slurm job for it is queued or running
#   STOPPED  CARP ran and exited without publishing: read the log named on the line
#   waiting  not started yet (the controller runs at most 4 at a time)
RP=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
REFS=/netscratch/dep_mercier/grp_marques/Aaryan/refs
cd "$RP" || exit 1
Q=$(squeue -u "$USER" -h -o '%i|%T|%M|%N|%k' 2>/dev/null)
for sp in $(awk -F, 'NR>1 && $1!="" {print $1}' config/annotation.csv); do
  man="$REFS/${sp}_repeats/MANIFEST.txt"; log="logs/annotate_repeats_carp/$sp.log"
  job=$(printf '%s\n' "$Q" | grep "wildcards_${sp}$" | head -1)
  if [ -s "$man" ]; then st="DONE      $(date -r "$man" '+%Y-%m-%d %H:%M')"
  elif [ -n "$job" ]; then st="RUNNING   $job"
  elif [ -s "$log" ] && grep -q "=== CARP exit" "$log"; then st="STOPPED   $(grep '=== CARP exit' "$log" | tail -1)  -> $RP/$log"
  else st="waiting"; fi
  printf '%-22s %s\n' "$sp" "$st"
done
L=$(ls -t logs/launch_carp_*.log 2>/dev/null | head -1)
if pgrep -u "$USER" -f "snakemake.*_repeats/MANIFEST.txt" > /dev/null; then c="running"; else c="finished or not running"; fi
echo "controller: $c   log: $RP/$L"
[ -n "$L" ] && tail -2 "$L"
