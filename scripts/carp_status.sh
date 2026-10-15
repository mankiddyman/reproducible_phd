#!/bin/bash
# Where the CARP runs are: Slurm state, CARP's steps (finished in any attempt, and in this one),
# last log lines, published or not, and for a running species an ETA (scripts/carp_eta.py).
# Usage: bash scripts/carp_status.sh            (from anywhere; a few seconds)
RP=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
REFS=/netscratch/dep_mercier/grp_marques/Aaryan/refs
PY=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk/bin/python
SPECIES="Drosera_paradoxa Drosera_binata"
echo "== Slurm jobs of $USER (rule_wildcards in the comment column)"
q=$(squeue -u "$USER" -o "%.10i %.9P %.45k %.8T %.11M %.5C %.8m %R" 2>/dev/null)
if echo "$q" | grep -q -E "repeats|carp"; then
    echo "$q" | grep -E "JOBID|repeats|carp" | sed 's/^/   /'
else
    echo "   no job with 'repeats' in its comment; all your jobs: $(echo "$q" | tail -n +2 | wc -l)"
fi
comments=$(squeue -u "$USER" -h -o "%k" 2>/dev/null)
for sp in $SPECIES; do
    log=$RP/logs/annotate_repeats_carp/$sp.log
    bm=$RP/results/$sp/repeats/carp/benchmarks
    echo "== $sp"
    if [ -s "$REFS/${sp}_repeats/MANIFEST.txt" ]; then
        echo "   PUBLISHED: $REFS/${sp}_repeats/ ($(stat -c %y "$REFS/${sp}_repeats/MANIFEST.txt" | cut -d. -f1))"
    fi
    n=$(ls "$bm"/*.tsv 2>/dev/null | wc -l)
    latest=$(ls -t "$bm"/*.tsv 2>/dev/null | head -n 3 | xargs -r -n 1 basename | sed 's/\.tsv$//' | paste -sd, -)
    echo "   CARP steps finished in any attempt: $n${latest:+ (latest: $latest)}"
    if [ -s "$log" ]; then
        head -n 2 "$log" | sed 's/^/   /'
        step=$(grep -E "of [0-9]+ steps \(" "$log" | tail -n 1)
        echo "   this attempt: ${step:-no step finished yet}"
        last=$(grep -E "^(rule|localrule) " "$log" | tail -n 1)
        echo "   last step started: ${last:-?}"
        grep -E "=== CARP exit|Error in rule|Exiting because" "$log" | tail -n 3 | sed 's/^/   /'
        echo "   log: $log"
    else
        echo "   no CARP log yet"
    fi
    # running now: an ETA scaled from a species whose CARP run has finished and been published
    if echo "$comments" | grep -q "wildcards_$sp" && [ -s "$log" ]; then
        for other in $SPECIES; do
            if [ "$other" != "$sp" ] && [ -s "$REFS/${other}_repeats/MANIFEST.txt" ]; then
                echo "   ETA, scaled from $other's finished run:"
                "$PY" "$RP/scripts/carp_eta.py" "$log" "$bm" "$RP/results/$other/repeats/carp/benchmarks" "${sp#Drosera_}" "${other#Drosera_}" 2>&1 \
                    | sed 's/^   /   /'
                break
            fi
        done
    fi
done
