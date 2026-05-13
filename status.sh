#!/bin/bash
# Status dashboard for the phd reproducible pipeline.
# Decodes snakemake's run-UUID job names by cross-referencing the driver log.

cd /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd

DRIVER_LOG=$(ls -t logs/sbatch/driver_*.out 2>/dev/null | head -1)

echo "=== $(date +%H:%M:%S) ==="
echo "Driver: $(squeue --me -h -o '%i %j %T %M' | grep phd_driver || echo NONE)"
echo

# Build a map: SLURM jobid -> (rule, species/wildcards)
declare -A JOB_RULE
declare -A JOB_WILD
if [ -f "$DRIVER_LOG" ]; then
    while IFS= read -r line; do
        # parse: "Job 24 has been submitted with SLURM jobid 788094 (log: .../rule_run_hifiasm/Drosera_filiformis/788094.log)"
        if [[ "$line" =~ SLURM\ jobid\ ([0-9]+)\ \(log:.*/rule_([^/]+)/([^/]+)/ ]]; then
            jobid="${BASH_REMATCH[1]}"
            rule="${BASH_REMATCH[2]}"
            wild="${BASH_REMATCH[3]}"
            JOB_RULE[$jobid]="$rule"
            JOB_WILD[$jobid]="$wild"
        fi
    done < "$DRIVER_LOG"
fi

echo "=== All my jobs (decoded) ==="
printf "%-10s %-8s %-10s %-30s %-25s %s\n" "JOBID" "STATE" "TIME" "RULE" "WILDCARD" "REASON"
while IFS=$'\t' read -r jobid state time reason; do
    rule="${JOB_RULE[$jobid]:-?}"
    wild="${JOB_WILD[$jobid]:-?}"
    printf "%-10s %-8s %-10s %-30s %-25s %s\n" "$jobid" "$state" "$time" "$rule" "$wild" "$reason"
done < <(squeue --me -h -o "%i	%T	%M	%R" | grep -v phd_driver)

echo
echo "=== Job stats by rule ==="
for jobid in "${!JOB_RULE[@]}"; do
    echo "${JOB_RULE[$jobid]}"
done | sort | uniq -c | sort -rn

echo
echo "=== Pipeline progress ==="
if [ -f "$DRIVER_LOG" ]; then
    grep -E "^[0-9]+ of [0-9]+ steps" "$DRIVER_LOG" | tail -1
fi

echo
echo "=== Disk ==="
df -h /netscratch/dep_mercier/grp_marques/Aaryan/ | tail -1 | awk '{print "  "$4" free, "$5" used"}'
