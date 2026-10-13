#!/bin/bash
# Where the CARP runs are: Slurm state, CARP's own step count, last log lines, published or not.
# Usage: bash scripts/carp_status.sh            (from anywhere; a few seconds)
RP=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
REFS=/netscratch/dep_mercier/grp_marques/Aaryan/refs
echo "== Slurm jobs of $USER (rule_wildcards in the comment column)"
q=$(squeue -u "$USER" -o "%.10i %.9P %.45k %.8T %.11M %.5C %.8m %R" 2>/dev/null)
if echo "$q" | grep -q -E "repeats|carp"; then
    echo "$q" | grep -E "JOBID|repeats|carp" | sed 's/^/   /'
else
    echo "   no job with 'repeats' in its comment; all your jobs: $(echo "$q" | tail -n +2 | wc -l)"
fi
for sp in Drosera_paradoxa Drosera_binata; do
    log=$RP/logs/annotate_repeats_carp/$sp.log
    echo "== $sp"
    if [ -s "$REFS/${sp}_repeats/MANIFEST.txt" ]; then
        echo "   PUBLISHED: $REFS/${sp}_repeats/ ($(stat -c %y "$REFS/${sp}_repeats/MANIFEST.txt" | cut -d. -f1))"
    fi
    if [ -s "$log" ]; then
        head -n 2 "$log" | sed 's/^/   /'
        step=$(grep -E "of [0-9]+ steps \(" "$log" | tail -n 1)
        echo "   CARP steps: ${step:-none finished yet}"
        last=$(grep -E "^(rule|localrule) " "$log" | tail -n 1)
        echo "   last rule started: ${last:-?}"
        grep -E "=== CARP exit|Error in rule|Exiting because" "$log" | tail -n 3 | sed 's/^/   /'
        echo "   log: $log"
    else
        echo "   no CARP log yet"
    fi
done
