#!/usr/bin/env bash
# carp_launch.sh -- CARP repeat + rDNA annotation (rule annotate_repeats_carp, then publish_repeats)
# for every species in config/annotation.csv not yet published to refs/. Checks first; submits only
# if every check passes:
#   0   this host resolves your user and Slurm accepts a submission from it (sbatch --test-only)
#   0b  every idle node in 'normal' can look up your user. A node that cannot (hpc07 on 2026-10-28:
#       "lookup userid 13760: connection refused") kills each CARP job in Singularity's first second,
#       is therefore always free, and attracts every job; it must be in the repeats exclude lists
#   1   Snakemake submits through Slurm (profiles/slurm); nothing else holds the repo
#   2   the CARP image is the one the published species were annotated with (sha256 in their MANIFESTs)
#   3   each job's Slurm time limit is >= 1.5x its runtime predicted from D. paradoxa (9.2 h for 4.30 Gb)
#   4   BAM/BAI/tabix in what CARP made for paradoxa (scorpioides has 2 chromosomes > 2^29 bp) -- informational
#   5   dry run restricted to the two CARP rules: published species get no jobs, every other species
#       exactly 2 (annotate + publish), no other rule, and every CARP job excludes the nodes that fail 0b.
#       The unrestricted dry run is reported alongside: under --rerun-triggers mtime, editing a config
#       file the assembly rules read (config/species.csv) makes Snakemake want to rebuild every frozen
#       assembly; --allowed-rules keeps CARP on the frozen FASTAs as they are.
# Then one detached controller (setsid nohup), at most 4 species at a time, scorpioides first.
# Report: logs/carp_gate_<time>.txt. Status at any time: scripts/carp_status.sh
RP=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
ENV=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk
SIF=/netscratch/dep_mercier/grp_marques/Aaryan/methods/carp/carp_1.9.5.sif
REFS=/netscratch/dep_mercier/grp_marques/Aaryan/refs
CFG=config/config.yaml
CARP_RULES="annotate_repeats_carp publish_repeats"
export PATH=$ENV/bin:$PATH
ME=$(id -u)
cd "$RP" || exit 1
TS=$(date +%Y%m%d_%H%M%S)
REP=$RP/logs/carp_gate_$TS.txt
STOPF=/tmp/carp_stop_$TS; : > "$STOPF"
stop() { echo "  STOP  $*"; echo "$*" >> "$STOPF"; }
SPECIES=$(awk -F, 'NR>1 && $1!="" {print $1}' config/annotation.csv)
jobs_tsv() {   # Snakemake dry-run log -> one line per planned job: rule <TAB> species <TAB> reason <TAB> Slurm --exclude
  awk 'function flush() { if (rule != "") print rule "\t" sp "\t" rsn "\t" ex; rule = ""; sp = ""; rsn = ""; ex = "" }
       /^(local)?(rule|checkpoint) [A-Za-z0-9_]+:/ { flush(); rule = $2; sub(/:$/, "", rule); next }
       /^[[:space:]]+reason:/ { rsn = $0; sub(/^[[:space:]]+reason:[[:space:]]*/, "", rsn); next }
       /^[[:space:]]+wildcards:/ { if (match($0, /species=[A-Za-z0-9_]+/)) sp = substr($0, RSTART + 8, RLENGTH - 8); next }
       /^[[:space:]]+resources:/ { if (match($0, /--exclude=[^ ]+/)) { ex = substr($0, RSTART + 10, RLENGTH - 10); sub(/,$/, "", ex) } next }
       /^Job stats:/ { flush() }
       END { flush() }' "$1"
}

{
echo "CARP on every species: checks first, submit only if all pass   ($(date '+%F %T'), $(hostname))"

echo; echo "=== 0. this host can run the controller"
if getent passwd "$ME" > /dev/null; then echo "  OK    $(hostname) resolves uid $ME"
else stop "$(hostname) cannot resolve your user (uid $ME). Run scripts/carp_launch.sh on hpc09"; fi
if sbatch --test-only -p normal -t 1 -c 1 --mem=100M --wrap true > /tmp/carp_sbtest_$TS.txt 2>&1; then
  echo "  OK    Slurm accepts a test submission from $(hostname)"
else stop "sbatch --test-only refused from $(hostname): $(head -1 /tmp/carp_sbtest_$TS.txt)"; fi

echo; echo "=== 0b. can each idle node in 'normal' look up your user?"
BAD=""; UNCHECKED=""; LINE=""
for n in $(sinfo -h -p normal -t idle,mix -o '%N' | xargs -r -n1 scontrol show hostnames 2> /dev/null | sort -u); do
  r=$(timeout 45 srun -I20 -Q -w "$n" -p normal -t 1 -c 1 --mem=100M -J carp_nodecheck \
        bash -c 'getent passwd "$(id -u)" > /dev/null && echo NSS_OK || echo NSS_FAIL' 2> /dev/null)
  case "$r" in
    *NSS_OK*)   LINE="$LINE $n:ok";;
    *NSS_FAIL*) LINE="$LINE $n:BAD"; BAD="$BAD $n";;
    *)          LINE="$LINE $n:?"; UNCHECKED="$UNCHECKED $n";;
  esac
done
echo "  $LINE" | fold -w 150 -s
if [ -n "$BAD" ]; then echo "  nodes that cannot look up uid $ME:$BAD (every CARP job must exclude them; checked in 5b)"
else echo "  OK    every idle node checked can look up your user"; fi
[ -n "$UNCHECKED" ] && echo "  note  not checked (no free slot within 20 s):$UNCHECKED"

echo; echo "=== 1. Snakemake submits to Slurm, nothing else holds the repo"
if grep -qE "^executor:[ \"']*slurm" profiles/slurm/config.yaml 2> /dev/null; then echo "  OK    profiles/slurm submits to Slurm"
else stop "profiles/slurm/config.yaml has no 'executor: slurm': jobs would run on $(hostname) itself"; fi
python3 -c "import snakemake_executor_plugin_slurm" 2> /dev/null && echo "  OK    slurm executor plugin, snakemake $(snakemake --version)" \
  || stop "snakemake-executor-plugin-slurm not importable from $ENV"
NL=$(ls .snakemake/locks 2> /dev/null | wc -l)
if [ "$NL" -gt 0 ]; then stop "$NL lock file(s) in .snakemake/locks: another Snakemake may be running in this repo ($(pgrep -u "$ME" -fc snakemake) snakemake processes on $(hostname); check the other nodes)"
else echo "  OK    no lock on the repo"; fi
if pgrep -u "$ME" -f "snakemake.*_repeats/MANIFEST.txt" > /dev/null; then stop "a CARP controller is already running on $(hostname): see scripts/carp_status.sh"; fi

echo; echo "=== 2. same image as the published species were annotated with"
if [ -s "$SIF" ]; then
  IMG=$(sha256sum "$SIF" | cut -d' ' -f1)
  REC=$(grep -hoE "[0-9a-f]{64}" "$REFS"/*_repeats/MANIFEST.txt 2> /dev/null | sort -u)
  if [ -z "$REC" ]; then echo "  note  no sha256 in the MANIFESTs; image not compared ($IMG)"
  elif echo "$REC" | grep -qx "$IMG"; then echo "  OK    sha256 $IMG matches the MANIFESTs"
  else stop "image sha256 $IMG is not the one in the MANIFESTs ($(echo $REC))"; fi
else stop "image missing: $SIF"; fi

echo; echo "=== 3. what each species asks Slurm for, against the runtime the finished runs predict"
python3 - "$CFG" "$REFS" << 'PY' || stop "a Slurm time limit is shorter than 1.5x the predicted runtime (table above)"
import csv, os, sys, yaml
CFG, REFS = sys.argv[1], sys.argv[2]
res = ((yaml.safe_load(open(CFG)) or {}).get("repeats") or {}).get("resources") or {}
RATE = 9.2 / 4.30   # h per Gb: D. paradoxa on 2026-10-15, the slower of the two finished runs
def hours(t):
    t, d = str(t), 0
    if "-" in t: d, t = t.split("-", 1); d = int(d)
    p = [int(x) for x in t.split(":")]
    if d: h = p[0] + (p[1] / 60 if len(p) > 1 else 0)
    elif len(p) == 3: h = p[0] + p[1] / 60
    else: h = p[0] / 60
    return 24 * d + h
bad = 0
print(f"  {'species':22s} {'state':9s} {'Gb':>5s} {'~h':>5s} {'limit h':>8s} {'thr':>4s} {'mem GB':>7s}")
for row in csv.DictReader(open("config/annotation.csv")):
    sp = row["species_id"]
    r = dict(res.get("default") or {}); r.update(res.get(sp) or {})
    fai = f"results/{sp}/assembly_final/{row['frozen_stage']}/{sp}_chr.fa.fai"
    gb = sum(int(l.split("\t")[1]) for l in open(fai)) / 1e9 if os.path.exists(fai) else float("nan")
    done = os.path.isfile(f"{REFS}/{sp}_repeats/MANIFEST.txt")
    try: lim = hours(r["time"])
    except Exception: lim = float("nan")
    ok = done or (lim >= 1.5 * RATE * gb)
    bad += not ok
    print(f"  {sp:22s} {'published' if done else 'to run':9s} {gb:5.2f} {RATE * gb:5.1f} {lim:8.1f} {str(r.get('threads', '?')):>4s} "
          f"{int(r.get('mem_mb', 0)) // 1000:7d}{'' if ok else '   <-- limit too short'}")
sys.exit(1 if bad else 0)
PY

echo; echo "=== 4. BAM/BAI/tabix in what CARP made for paradoxa (scorpioides: 2 chromosomes > 2^29 bp)"
STEPS=$(grep -ohE "^(local)?rule [A-Za-z0-9_]+" logs/annotate_repeats_carp/Drosera_paradoxa.log 2> /dev/null | awk '{print $2}' | sort -u)
IDX=$(find results/Drosera_paradoxa/repeats -type f \( -name "*.bai" -o -name "*.csi" -o -name "*.tbi" -o -name "*.bam" \) 2> /dev/null | head -5 | tr '\n' ' ')
ISTEP=$(echo "$STEPS" | grep -iE "bam|bai|tabix|index" | tr '\n' ' ')
if [ -z "$STEPS" ]; then echo "  note  no step list in the paradoxa log, so unverified"
elif [ -z "$IDX$ISTEP" ]; then echo "  OK    none of CARP's $(echo $STEPS | wc -w) steps or files is a BAM/BAI/tabix index: the 2^29 limit does not touch scorpioides"
else echo "  note  CARP makes indexes ($ISTEP$IDX); a failure there is a side step, which the rule tolerates once the unified GFF and its statistics exist"; fi

ALL=$(for sp in $SPECIES; do printf '%s ' "$REFS/${sp}_repeats/MANIFEST.txt"; done)
echo; echo "=== 5a. what an unrestricted run would also redo (reported only)"
snakemake $ALL -n --nolock --rerun-triggers mtime --profile profiles/slurm > logs/carp_dry_full_$TS.log 2>&1
jobs_tsv logs/carp_dry_full_$TS.log > /tmp/carp_jobs_full_$TS.tsv
echo "  $(wc -l < /tmp/carp_jobs_full_$TS.tsv) jobs, of which outside the CARP rules: $(awk -F'\t' '$1 != "annotate_repeats_carp" && $1 != "publish_repeats"' /tmp/carp_jobs_full_$TS.tsv | wc -l)"
echo "  root causes outside the CARP rules:"
awk -F'\t' '$1 != "annotate_repeats_carp" && $1 != "publish_repeats" && $3 !~ /^Input files updated by another job/ {print $1 ": " substr($3, 1, 150)}' \
  /tmp/carp_jobs_full_$TS.tsv | sort | uniq -c | sort -rn | head -8 | sed 's/^/   /'

echo; echo "=== 5b. dry run restricted to the CARP rules (the check)"
snakemake $ALL --allowed-rules $CARP_RULES -n --nolock --rerun-triggers mtime --profile profiles/slurm > logs/carp_dry_$TS.log 2>&1
rc=$?
awk '/^Job stats:/{f=1} f{print} f && /^total/{exit}' logs/carp_dry_$TS.log | sed 's/^/  /'
if [ $rc -ne 0 ]; then
  stop "restricted dry run failed (exit $rc), log $RP/logs/carp_dry_$TS.log"
  grep -m1 -A4 -E "Exception|Error" logs/carp_dry_$TS.log | sed 's/^/    /'
else
  jobs_tsv logs/carp_dry_$TS.log > /tmp/carp_jobs_$TS.tsv
  OTHER=$(awk -F'\t' '$1 != "annotate_repeats_carp" && $1 != "publish_repeats" {print $1}' /tmp/carp_jobs_$TS.tsv | sort -u | tr '\n' ' ')
  [ -n "$OTHER" ] && stop "rules other than CARP's would run: $OTHER"
  for sp in $SPECIES; do
    n=$(awk -F'\t' -v s="$sp" '$2 == s' /tmp/carp_jobs_$TS.tsv | wc -l)
    ex=$(awk -F'\t' -v s="$sp" '$1 == "annotate_repeats_carp" && $2 == s {print $4; exit}' /tmp/carp_jobs_$TS.tsv)
    if [ -s "$REFS/${sp}_repeats/MANIFEST.txt" ]; then
      if [ "$n" -eq 0 ]; then echo "  OK    $sp: published, no jobs"
      else stop "$sp is published but would rerun: $(awk -F'\t' -v s="$sp" '$2 == s {print $1 " (" substr($3, 1, 120) ")"; exit}' /tmp/carp_jobs_$TS.tsv)"; fi
      continue
    fi
    if [ "$n" -ne 2 ]; then stop "$sp: restricted dry run plans $n jobs, expected 2"; continue; fi
    miss=""; for b in $BAD; do case ",$ex," in *",$b,"*) ;; *) miss="$miss $b";; esac; done
    if [ -n "$miss" ]; then stop "$sp's CARP job would not exclude$miss (its exclude: ${ex:-none})"
    else echo "  OK    $sp: annotate + publish, excluding ${ex:-nothing}"; fi
  done
fi
} 2>&1 | tee "$REP"

NEW=""; for sp in $SPECIES; do [ -s "$REFS/${sp}_repeats/MANIFEST.txt" ] || NEW="$NEW $sp"; done
if [ -s "$STOPF" ]; then
  echo; echo "NOT LAUNCHED: $(wc -l < "$STOPF") check(s) failed:"; sed 's/^/  - /' "$STOPF"
  echo "Paste back $REP"
elif [ -z "$NEW" ]; then
  echo; echo "Nothing to launch: every species is published."
else
  LOG=logs/launch_carp_$TS.log
  NT=$(for sp in $NEW; do printf '%s ' "$REFS/${sp}_repeats/MANIFEST.txt"; done)
  PRIO=""
  case " $NEW " in *" Drosera_scorpioides "*)
    snakemake --help 2> /dev/null | grep -q -- "--prioritize" && PRIO="--prioritize $REFS/Drosera_scorpioides_repeats/MANIFEST.txt";;
  esac
  setsid nohup bash -c "snakemake $NT --allowed-rules $CARP_RULES --profile profiles/slurm --rerun-triggers mtime --keep-going --jobs 4 $PRIO; echo \"snakemake exit \$? at \$(date)\"" > "$LOG" 2>&1 < /dev/null &
  echo; echo "LAUNCHED:$NEW"
  echo "controller log: $RP/$LOG   checks: $REP"
  sleep 90
  "$RP/scripts/carp_status.sh"
fi
