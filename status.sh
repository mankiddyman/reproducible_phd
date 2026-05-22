#!/bin/bash
# PHD pipeline status dashboard.
# Shows per-species progression through hifiasm → blobtoolkit → decontamination
# → QC, with ploidy-aware target selection. Surfaces failures and bottlenecks.
set -uo pipefail
cd /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd

DRIVER_LOG=$(ls -t logs/sbatch/driver_*.out 2>/dev/null | head -1)

# ---- Colors ----
G="\033[32m"; R="\033[31m"; Y="\033[33m"; B="\033[34m"; C="\033[36m"
D="\033[2m"; W="\033[1m"; N="\033[0m"

echo
echo -e "${W}═════════════════════════════════════════════════════════════════${N}"
echo -e "${W}  PHD pipeline status  —  $(date '+%a %H:%M:%S')${N}"
echo -e "${W}═════════════════════════════════════════════════════════════════${N}"

# ---- Driver ----
driver_line=$(squeue --me -h -o "%i %j %T %M" | awk '$2 == "phd_driver"')
if [ -n "$driver_line" ]; then
    read -r drv_id _ drv_state drv_time <<< "$driver_line"
    echo -e "Driver: ${G}${drv_id}${N} ${drv_state} (running ${drv_time})"
else
    drv_id=""
    echo -e "Driver: ${R}NONE${N}  ${D}(restart: sbatch scripts/launchers/launch_phd_redo.sh)${N}"
fi

# ---- Build job rule/wildcards lookup ----
declare -A JOB_RULE JOB_WILD
if [ -d snakemake_tmp/slurm_logs ]; then
    while IFS= read -r path; do
        if [[ "$path" =~ slurm_logs/rule_([^/]+)/([^/]+)/([0-9]+)\.log$ ]]; then
            JOB_RULE["${BASH_REMATCH[3]}"]="${BASH_REMATCH[1]}"
            JOB_WILD["${BASH_REMATCH[3]}"]="${BASH_REMATCH[2]}"
        fi
    done < <(find snakemake_tmp/slurm_logs -name "*.log" -mtime -3 2>/dev/null)
fi

# ---- Active jobs (compact) ----
echo
echo -e "${W}── Active jobs${N}"
mapfile -t jobs < <(squeue --me -h -o "%i|%T|%M|%R")
if [ ${#jobs[@]} -eq 0 ]; then
    echo -e "  ${D}(none — driver is idle or stopped)${N}"
else
    printf "  %-9s %-9s %-9s %-10s %-32s %s\n" "JOBID" "STATE" "TIME" "NODE" "RULE" "WILDCARDS"
    declare -A RULE_COUNT
    for line in "${jobs[@]}"; do
        IFS='|' read -r id state time node <<< "$line"
        rule="${JOB_RULE[$id]:-?}"
        wild="${JOB_WILD[$id]:-?}"
        if [ "$id" = "$drv_id" ]; then rule="(driver)"; wild="-"; fi
        # Detect nextflow children by name pattern from squeue
        if [ "$rule" = "?" ]; then
            jname=$(squeue -j "$id" -h -o "%j" 2>/dev/null)
            if [[ "$jname" =~ ^nf- ]]; then
                rule="nextflow"
                wild=$(echo "$jname" | sed 's/^nf-SANGERTOL_BLOBTOOLKIT_//' | head -c 40)
            fi
        fi
        # Highlight long-running based on rule
        time_color=""
        # Parse time as minutes (handles HH:MM:SS or D-HH:MM:SS)
        tmin=$(echo "$time" | awk -F'[-:]' '{
            if (NF==4) { print $1*1440 + $2*60 + $3 + $4/60 }
            else if (NF==3) { print $1*60 + $2 + $3/60 }
            else if (NF==2) { print $1 + $2/60 }
            else { print 0 }
        }')
        # Warn if hifiasm > 24h or blobtoolkit > 48h or anything > 7d
        case "$rule" in
            run_hifiasm) if (( $(echo "$tmin > 1440" | bc -l) )); then time_color="$Y"; fi ;;
            run_blobtoolkit_initial|nextflow) if (( $(echo "$tmin > 2880" | bc -l) )); then time_color="$Y"; fi ;;
        esac
        if (( $(echo "$tmin > 10080" | bc -l) )); then time_color="$R"; fi
        case "$state" in
            RUNNING) sc="${G}RUN${N}" ;;
            PENDING) sc="${Y}PEND${N}" ;;
            *) sc="${state:0:8}" ;;
        esac
        printf "  %-9s %-15b %-9b %-10s %-32s %s\n" "$id" "$sc" "${time_color}${time}${N}" "$node" "$rule" "$wild"
        if [ "$rule" != "?" ] && [ "$rule" != "(driver)" ]; then
            RULE_COUNT["$rule"]=$((${RULE_COUNT["$rule"]:-0}+1))
        fi
    done

    if [ ${#RULE_COUNT[@]} -gt 0 ]; then
        echo
        echo -e "  ${D}counts:${N} $(for k in "${!RULE_COUNT[@]}"; do printf "%s=%d  " "$k" "${RULE_COUNT[$k]}"; done)"
    fi
fi

# ---- Per-species stage tracker (ploidy-aware) ----
echo
echo -e "${W}── Per-species pipeline progress${N}"

# Parse ploidy from species.csv to determine relevant targets
declare -A PLOIDY
if [ -f config/species.csv ]; then
    awk -F, 'NR>1 {gsub(/"/,""); print $1"="$NF}' config/species.csv 2>/dev/null | while read -r kv; do
        echo "$kv"
    done > /tmp/.ploidy_$$
fi

# Use a python helper for ploidy lookup since csv parsing in bash is fragile
get_ploidy() {
    python3 -c "
import csv, sys
with open('config/species.csv') as f:
    for row in csv.DictReader(f):
        if row['species_id'] == '$1':
            v = row.get('exp_ploidy', '').strip()
            try: print(int(v))
            except: print(2)
            break
    else:
        print(2)
"
}

# Stage indicators
stage_status() {
    # Args: species, target (hap1/hap2/p_utg)
    # Returns colored emoji-style summary: H S F B D Q
    local sp=$1 tgt=$2
    local hifiasm_done="results/${sp}/hifiasm/${sp}.done"
    local manifest="results/${sp}/assembly/initial/manifest.tsv"
    local fasta="results/${sp}/assembly/initial/${tgt}/${sp}.fa"
    local blob_done="results/${sp}/blobtoolkit/initial/${tgt}/.done"
    local blob_params="results/${sp}/blobtoolkit/initial/${tgt}/params.yaml"
    local decontam="results/${sp}/assembly/decontaminated/${tgt}/${sp}.fa"
    local qc_compleasm="results/${sp}/qc/compleasm/decontaminated/${tgt}/summary.txt"

    local cells=()
    # H — hifiasm
    if [ -f "$hifiasm_done" ]; then cells+=("${G}H${N}"); else cells+=("${D}h${N}"); fi
    # S — standardize
    if [ -f "$manifest" ]; then cells+=("${G}S${N}"); else cells+=("${D}s${N}"); fi
    # F — fasta available for this target
    if [ -f "$fasta" ] && [ -s "$fasta" ]; then cells+=("${G}F${N}"); else cells+=("${D}f${N}"); fi
    # B — blobtoolkit
    if [ -f "$blob_done" ]; then
        cells+=("${G}B${N}")
    elif [ -f "$blob_params" ]; then
        cells+=("${Y}B${N}")
    else
        cells+=("${D}b${N}")
    fi
    # D — decontamination
    if [ -f "$decontam" ]; then cells+=("${G}D${N}"); else cells+=("${D}d${N}"); fi
    # Q — qc on decontaminated
    if [ -f "$qc_compleasm" ]; then cells+=("${G}Q${N}"); else cells+=("${D}q${N}"); fi

    echo "${cells[*]}" | tr -d ' '
}

# Legend
echo -e "  ${D}H${N}=hifiasm  ${D}S${N}=standardize  ${D}F${N}=fasta  ${D}B${N}=blobtoolkit  ${D}D${N}=decontam  ${D}Q${N}=qc"
echo -e "  ${D}${G}green${N}=done  ${D}${Y}yellow${N}=in-progress  ${D}${D}grey${N}=pending${N}"
echo
printf "  %-22s %-8s %-12s %-12s\n" "SPECIES" "PLOIDY" "TARGETS" "PROGRESS"

for sp_dir in results/Drosera_*; do
    [ -d "$sp_dir" ] || continue
    sp=$(basename "$sp_dir")

    ploidy=$(get_ploidy "$sp")
    if [ "$ploidy" = "2" ] || [ -z "$ploidy" ]; then
        targets=("hap1" "hap2")
    else
        targets=("p_utg")
    fi

    target_str=$(IFS=','; echo "${targets[*]}")
    progress_strs=()
    for tgt in "${targets[@]}"; do
        progress_strs+=("$tgt:$(stage_status "$sp" "$tgt")")
    done

    printf "  %-22s %-8s %-12s " "$sp" "$ploidy" "$target_str"
    for ps in "${progress_strs[@]}"; do
        echo -ne "$ps  "
    done
    echo
done

# ---- Long-running jobs warning ----
echo
echo -e "${W}── Health checks${N}"
warnings=()

# Check for hifiasm > 24h
for line in "${jobs[@]:-}"; do
    IFS='|' read -r id state time _ <<< "$line"
    rule="${JOB_RULE[$id]:-}"
    if [ "$rule" = "run_hifiasm" ]; then
        tmin=$(echo "$time" | awk -F'[-:]' '{
            if (NF==4) { print $1*1440 + $2*60 + $3 + $4/60 }
            else if (NF==3) { print $1*60 + $2 + $3/60 }
            else { print 0 }
        }')
        if (( $(echo "$tmin > 1440" | bc -l) )); then
            warnings+=("${Y}⚠${N}  hifiasm job $id (${JOB_WILD[$id]:-?}) running >24h — check log")
        fi
    fi
done

# Recent OOMs from sacct
oom_count=$(sacct -X --starttime now-1day --format=State -P 2>/dev/null | grep -c "OUT_OF_MEMORY" || true)
oom_count=$(echo "${oom_count:-0}" | tr -d '\n' | awk '{print $1+0}')
if [ "$oom_count" -gt 0 ]; then
    warnings+=("${Y}⚠${N}  $oom_count OOM event(s) in past 24h (sacct --state=OUT_OF_MEMORY for details)")
fi

# Recent failures from driver log
if [ -f "$DRIVER_LOG" ]; then
    recent_fails=$(grep -E "SLURM-job '[0-9]+' failed" "$DRIVER_LOG" 2>/dev/null | tail -5)
    if [ -n "$recent_fails" ]; then
        n_fails=$(echo "$recent_fails" | wc -l)
        warnings+=("${Y}⚠${N}  $n_fails recent failure(s) in driver log")
        while IFS= read -r fline; do
            fid=$(echo "$fline" | grep -oE "'[0-9]+'" | tr -d "'")
            frule="${JOB_RULE[$fid]:-?}"
            fwild="${JOB_WILD[$fid]:-?}"
            warnings+=("    ${R}✗${N} $fid  $frule  $fwild")
        done <<< "$recent_fails"
    fi
fi

# Driver missing while jobs running
if [ -z "$drv_id" ] && [ ${#jobs[@]} -gt 0 ]; then
    warnings+=("${R}!${N}  Driver DEAD but ${#jobs[@]} orphan child jobs still running — they will complete but no new jobs will submit")
fi

# Driver alive but only 3 submissions (known throttle bug)
if [ -n "$drv_id" ] && [ -n "$DRIVER_LOG" ] && [ -f "$DRIVER_LOG" ]; then
    n_submit=$(grep -c "has been submitted with SLURM jobid" "$DRIVER_LOG" 2>/dev/null || echo 0)
    if [ "$n_submit" -le 4 ]; then
        # Get driver age in minutes
        drv_min=$(echo "$drv_time" | awk -F'[-:]' '{
            if (NF==4) { print $1*1440 + $2*60 + $3 + $4/60 }
            else if (NF==3) { print $1*60 + $2 + $3/60 }
            else { print $1 + $2/60 }
        }')
        if (( $(echo "$drv_min > 30" | bc -l) )); then
            warnings+=("${Y}⚠${N}  Driver only submitted $n_submit jobs in ${drv_time} — scheduler-throttle bug likely")
        fi
    fi
fi

if [ ${#warnings[@]} -eq 0 ]; then
    echo -e "  ${G}✓${N}  no issues detected"
else
    for w in "${warnings[@]}"; do
        echo -e "  $w"
    done
fi

# ---- Snakemake progress ----
if [ -f "$DRIVER_LOG" ]; then
    last_pct=$(grep -E "^[0-9]+ of [0-9]+ steps" "$DRIVER_LOG" 2>/dev/null | tail -1)
    if [ -n "$last_pct" ]; then
        n_submit=$(grep -c "has been submitted with SLURM jobid" "$DRIVER_LOG" 2>/dev/null || echo 0)
        echo
        echo -e "${W}── Pipeline progress${N}"
        echo "  $last_pct"
        echo "  driver submissions so far: $n_submit"
    fi
fi

# ---- Disk ----
echo
df -h /netscratch/dep_mercier/grp_marques 2>/dev/null | tail -1 | awk '{
    pct=$5; gsub("%","",pct);
    if (pct+0 > 95) color="\033[31m"; else if (pct+0 > 90) color="\033[33m"; else color="\033[32m";
    printf "Disk: %s%s%%\033[0m used  (%s free of %s)\n", color, pct, $4, $2
}'

# Cleanup
rm -f /tmp/.ploidy_$$
