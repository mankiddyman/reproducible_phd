#!/bin/bash
# Status dashboard for the phd reproducible pipeline.
# Decodes SLURM job names by finding their log path in snakemake_tmp/slurm_logs/
# (more reliable than sacct --comment, which our cluster doesn't expose).
set -uo pipefail
cd /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd

DRIVER_LOG=$(ls -t logs/sbatch/driver_*.out 2>/dev/null | head -1)

# ANSI colors
G="\033[32m"; R="\033[31m"; Y="\033[33m"; D="\033[2m"; N="\033[0m"

echo "=========================================="
echo " PHD pipeline status  —  $(date +%H:%M:%S)"
echo "=========================================="

# --- driver ---
driver=$(squeue --me -h -o "%i %j %T %M" | awk '$2 == "phd_driver"')
if [ -n "$driver" ]; then
    echo -e "Driver: ${G}${driver}${N}"
else
    echo -e "Driver: ${R}NONE${N}"
fi

# --- active jobs (decoded via slurm_logs path) ---
echo
echo "--- Active jobs ---"
declare -A JOB_RULE JOB_WILD

# Build lookup from snakemake_tmp/slurm_logs/
if [ -d snakemake_tmp/slurm_logs ]; then
    while IFS= read -r path; do
        # Path: snakemake_tmp/slurm_logs/rule_<rule>/<wildcards>/<jobid>.log
        if [[ "$path" =~ slurm_logs/rule_([^/]+)/([^/]+)/([0-9]+)\.log$ ]]; then
            JOB_RULE["${BASH_REMATCH[3]}"]="${BASH_REMATCH[1]}"
            JOB_WILD["${BASH_REMATCH[3]}"]="${BASH_REMATCH[2]}"
        fi
    done < <(find snakemake_tmp/slurm_logs -name "*.log" -mtime -2 2>/dev/null)
fi

mapfile -t jobs < <(squeue --me -h -o "%i|%T|%M|%R")
if [ ${#jobs[@]} -eq 0 ]; then
    echo "  (none)"
else
    printf "%-10s %-17s %-10s %-10s %-32s %s\n" "JOBID" "STATE" "TIME" "NODE" "RULE" "WILDCARDS"
    declare -A RULE_COUNT
    for line in "${jobs[@]}"; do
        IFS='|' read -r id state time node <<< "$line"
        rule="${JOB_RULE[$id]:-?}"
        wild="${JOB_WILD[$id]:-?}"
        # Filter known special jobs
        if [ "$id" = "${driver%% *}" ]; then
            rule="(driver)"; wild="-"
        fi
        # Color state
        case "$state" in
            RUNNING) sc="${G}${state}${N}" ;;
            PENDING) sc="${Y}${state}${N}" ;;
            *) sc="${state}" ;;
        esac
        printf "%-10s ${sc}%-$((17-${#state}))s${N} %-10s %-10s %-32s %s\n" \
            "$id" "" "$time" "$node" "$rule" "$wild"
        if [ "$rule" != "?" ] && [ "$rule" != "(driver)" ]; then
            RULE_COUNT["$rule"]=$((${RULE_COUNT["$rule"]:-0}+1))
        fi
    done

    # --- counts by rule ---
    echo
    echo "--- Active rules ---"
    if [ ${#RULE_COUNT[@]} -eq 0 ]; then
        echo "  (none decoded)"
    else
        for k in "${!RULE_COUNT[@]}"; do
            printf "  %-32s %d\n" "$k" "${RULE_COUNT[$k]}"
        done | sort
    fi
fi

# --- per-species pipeline state ---
echo
echo "--- Per-species pipeline state ---"
printf "%-26s %-10s %-12s %s\n" "SPECIES" "HIFIASM" "STANDARDIZE" "BLOBTOOLKIT"
for sp_dir in results/Drosera_*; do
    [ -d "$sp_dir" ] || continue
    sp=$(basename "$sp_dir")

    # Hifiasm done?
    if [ -f "$sp_dir/hifiasm/$sp.done" ]; then
        hif="${G}done${N}"
    else
        hif="${R}--${N}"
    fi

    # Standardize done?
    if [ -f "$sp_dir/assembly/initial/manifest.tsv" ]; then
        std="${G}done${N}"
    else
        std="${R}--${N}"
    fi

    # Blobtoolkit state per target
    blob_parts=()
    for blob_dir in "$sp_dir"/blobtoolkit/initial/*/; do
        [ -d "$blob_dir" ] || continue
        target=$(basename "$blob_dir")
        if [ -f "$blob_dir/.done" ]; then
            blob_parts+=("${G}${target}✓${N}")
        elif [ -f "$blob_dir/params.yaml" ]; then
            blob_parts+=("${Y}${target}…${N}")
        else
            blob_parts+=("${D}${target}-${N}")
        fi
    done
    blob_str=$(IFS=' '; echo "${blob_parts[*]:-${D}none${N}}")
    printf "%-26s ${hif}%-$((10-4))s${N} ${std}%-$((12-${#std}+4))s${N} %b\n" "$sp" "" "" "$blob_str"
done

# --- recent failures from driver log ---
echo
echo "--- Recent failures (last 5) ---"
if [ -f "$DRIVER_LOG" ]; then
    fails=$(grep -E "SLURM-job '[0-9]+' failed" "$DRIVER_LOG" 2>/dev/null \
        | sed -E "s/.*SLURM-job '([0-9]+)' failed.*/\1/" | sort -u | tail -5)
    if [ -n "$fails" ]; then
        while read -r fail_id; do
            rule="${JOB_RULE[$fail_id]:-?}"
            wild="${JOB_WILD[$fail_id]:-?}"
            echo "  $fail_id  $rule  $wild"
        done <<< "$fails"
    else
        echo "  (none)"
    fi
fi

# --- pipeline progress ---
echo
echo "--- Pipeline progress ---"
if [ -f "$DRIVER_LOG" ]; then
    grep -E "^[0-9]+ of [0-9]+ steps" "$DRIVER_LOG" 2>/dev/null | tail -1 | sed 's/^/  /'
fi

# --- disk ---
echo
df -h /netscratch/dep_mercier/grp_marques 2>/dev/null | tail -1 | awk '{printf "Disk: %s free / %s used (%s)\n", $4, $3, $5}'
