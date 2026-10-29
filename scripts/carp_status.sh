#!/usr/bin/env bash
# carp_status.sh -- where every species' CARP repeat annotation stands, one line each.
#   DONE        published (MANIFEST.txt written by rule publish_repeats)
#   RUNNING / PENDING   its Slurm job: job id, elapsed time, node
#   ANNOTATED   CARP finished (carp.done); publish_repeats still to run
#   FAILED      its job ended in the current launch without publishing: read the log named on the line
#   waiting     not started in the current launch (the controller runs at most 4 at a time)
# "Current launch" = the newest logs/launch_carp_*.log; log lines from earlier launches are ignored.
RP=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
REFS=/netscratch/dep_mercier/grp_marques/Aaryan/refs
ME=$(id -u)
cd "$RP" || exit 1
L=$(ls -t logs/launch_carp_*.log 2> /dev/null | head -1)
T0=0
if [ -n "$L" ]; then
  s=$(basename "$L" .log); s=${s#launch_carp_}
  hh=${s:9:2}; mm=${s:11:2}; ss=${s:13:2}; [ -z "$ss" ] && ss=00
  T0=$(date -d "${s:0:8} $hh:$mm:$ss" +%s 2> /dev/null || echo 0)
fi
Q=$(squeue -u "$ME" -h -o '%i|%T|%M|%N|%k' 2> /dev/null)
for sp in $(awk -F, 'NR>1 && $1!="" {print $1}' config/annotation.csv); do
  man="$REFS/${sp}_repeats/MANIFEST.txt"; log="logs/annotate_repeats_carp/$sp.log"
  done_f="results/$sp/repeats/carp_run/carp.done"
  job=$(printf '%s\n' "$Q" | grep "wildcards_${sp}$" | head -1)
  if [ -s "$man" ]; then st="DONE        published $(date -r "$man" '+%Y-%m-%d %H:%M')"
  elif [ -n "$job" ]; then
    IFS='|' read -r jid jst jel jnode jc <<< "$job"
    st=$(printf '%-11s job %s, %s, %s' "$jst" "$jid" "$jel" "${jnode:-no node yet}")
    case "$jc" in rule_publish_repeats*) st="$st (publishing)";; esac
  elif [ -s "$done_f" ]; then st="ANNOTATED   carp.done written; publish_repeats still to run"
  elif [ -s "$log" ] && [ "$(date -r "$log" +%s)" -ge "$T0" ]; then
    st="FAILED      $(grep '=== CARP exit' "$log" | tail -1 | sed 's/^=== //')  -> $RP/$log"
  else st="waiting"; fi
  printf '%-22s %s\n' "$sp" "$st"
done
if pgrep -u "$ME" -f "carp_controller_|snakemake.*_repeats/MANIFEST.txt" > /dev/null; then c="running on $(hostname)"
else c="not running on $(hostname) (if you launched elsewhere, check there)"; fi
echo "controller: $c   log: $RP/$L"
[ -n "$L" ] && tail -2 "$L"
