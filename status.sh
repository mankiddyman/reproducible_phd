#!/bin/bash
echo "=== $(date +%H:%M:%S) ==="
echo "Driver: $(squeue --me -h -o '%j %T %M' | grep phd_ || echo NONE)"
echo "Children: $(squeue --me -h | grep nf-SAN | wc -l) jobs (R + PD)"
echo
echo "Per-species (latest active step):"
for sp in Drosera_binata Drosera_paradoxa Drosera_roseana Drosera_scorpioides Drosera_aliciae Drosera_tokaiensis Drosera_filiformis; do
  log="logs/blobtoolkit/${sp}_run.log"
  running=$(squeue --me -h -o "%j" | grep -c "${sp}_" 2>/dev/null)
  done_marker=""
  [ -f "results/${sp}/blobtoolkit/initial/.done" ] && done_marker=" [DONE]"
  # Get from SLURM queue what's actually running for this species
  # Strip the long prefix to get the meaningful tail of the process name
  current=$(squeue --me -h -o "%j" | grep "${sp}_" | head -1 | sed -E 's/.*BLOBTOOLKIT_//; s/_\(.*//')
  if [ -z "$current" ]; then
    # No active SLURM job — look at last line of log
    current=$(tail -200 "$log" 2>/dev/null | grep -E "^\[[a-f0-9_-]+/" | tail -1 | sed -E 's/^\[[^]]+\] //; s/^SAN[^:]+://; s/ \| .*//' | tr -d '…')
  fi
  printf "  %-22s active=%d  stage=%s%s\n" "${sp#Drosera_}" "$running" "${current:-idle}" "$done_marker"
done
echo
echo "Disk: $(df -h /netscratch/dep_mercier/grp_marques/Aaryan/ | tail -1 | awk '{print $4 " free, " $5 " used"}')"
echo "Failures total: $(find /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/nextflow_work -name .exitcode 2>/dev/null -exec sh -c 'ec=$(cat "$1" 2>/dev/null); [ "$ec" != "0" ] && [ -n "$ec" ] && echo y' _ {} \; | wc -l)"
echo "Failures last hour: $(find /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/nextflow_work -name .exitcode -mmin -60 2>/dev/null -exec sh -c 'ec=$(cat "$1" 2>/dev/null); [ "$ec" != "0" ] && [ -n "$ec" ] && echo y' _ {} \; | wc -l)"
