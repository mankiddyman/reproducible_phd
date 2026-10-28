#!/usr/bin/env bash
# carp_launch.sh -- CARP repeat + rDNA annotation (rule annotate_repeats_carp, then publish_repeats)
# for every species in config/annotation.csv not yet published to refs/. Checks first; submits only
# if every check passes:
#   1  Snakemake submits through Slurm (profiles/slurm); nothing else holds the repo
#   2  the CARP image is the one the published species were annotated with (sha256 in their MANIFESTs)
#   3  each job's Slurm time limit is >= 1.5x its runtime predicted from D. paradoxa (9.2 h for 4.30 Gb)
#   4  BAM/BAI/tabix in what CARP made for paradoxa (scorpioides has 2 chromosomes > 2^29 bp,
#      beyond what BAI or .tbi can address) -- informational
#   5  dry run: published species get no jobs, every other species exactly 2 (annotate + publish)
# Then one detached controller (setsid nohup): at most 4 species at a time, scorpioides first.
# Report: /tmp/carp_gate_<time>.txt. Status at any time: scripts/carp_status.sh
RP=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
ENV=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk
SIF=/netscratch/dep_mercier/grp_marques/Aaryan/methods/carp/carp_1.9.5.sif
REFS=/netscratch/dep_mercier/grp_marques/Aaryan/refs
CFG=config/config.yaml
export PATH=$ENV/bin:$PATH
cd "$RP" || exit 1
TS=$(date +%Y%m%d_%H%M%S)
REP=/tmp/carp_gate_$TS.txt
STOPF=/tmp/carp_stop_$TS; : > "$STOPF"
stop() { echo "  STOP  $*"; echo "$*" >> "$STOPF"; }
SPECIES=$(awk -F, 'NR>1 && $1!="" {print $1}' config/annotation.csv)

{
echo "CARP on every species: checks first, submit only if all pass   ($(date '+%F %T'), $(hostname))"

echo; echo "=== 1. Snakemake submits to Slurm, nothing else holds the repo"
if grep -qE "^executor:[ \"']*slurm" profiles/slurm/config.yaml 2>/dev/null; then echo "  OK    profiles/slurm submits to Slurm"
else stop "profiles/slurm/config.yaml has no 'executor: slurm': jobs would run on $(hostname) itself"; fi
python3 -c "import snakemake_executor_plugin_slurm" 2>/dev/null && echo "  OK    slurm executor plugin, snakemake $(snakemake --version)" \
  || stop "snakemake-executor-plugin-slurm not importable from $ENV"
NL=$(ls .snakemake/locks 2>/dev/null | wc -l)
if [ "$NL" -gt 0 ]; then stop "$NL lock file(s) in .snakemake/locks: another Snakemake may be running in this repo ($(pgrep -u "$USER" -fc snakemake) snakemake processes on $(hostname); check the other nodes)"
else echo "  OK    no lock on the repo"; fi
if pgrep -u "$USER" -f "snakemake.*_repeats/MANIFEST.txt" > /dev/null; then stop "a CARP controller is already running on $(hostname): see scripts/carp_status.sh"; fi

echo; echo "=== 2. same image as the published species were annotated with"
if [ -s "$SIF" ]; then
  IMG=$(sha256sum "$SIF" | cut -d' ' -f1)
  REC=$(grep -hoE "[0-9a-f]{64}" "$REFS"/*_repeats/MANIFEST.txt 2>/dev/null | sort -u)
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
print(f"  {'species':22s} {'state':9s} {'Gb':>5s} {'~h':>5s} {'limit h':>8s} {'thr':>4s} {'mem GB':>7s}  exclude")
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
          f"{int(r.get('mem_mb', 0)) // 1000:7d}  {r.get('exclude', 'hpc15')}{'' if ok else '   <-- limit too short'}")
sys.exit(1 if bad else 0)
PY

echo; echo "=== 4. BAM/BAI/tabix in what CARP made for paradoxa (scorpioides: 2 chromosomes > 2^29 bp)"
STEPS=$(grep -ohE "^(local)?rule [A-Za-z0-9_]+" logs/annotate_repeats_carp/Drosera_paradoxa.log 2>/dev/null | awk '{print $2}' | sort -u)
echo "  CARP steps that ran ($(echo $STEPS | wc -w)): $(echo $STEPS | tr ' ' ',')" | fold -w 150 -s | sed '2,$s/^/    /'
IDX=$(find results/Drosera_paradoxa/repeats -type f \( -name "*.bai" -o -name "*.csi" -o -name "*.tbi" -o -name "*.bam" \) 2>/dev/null | head -5 | tr '\n' ' ')
ISTEP=$(echo "$STEPS" | grep -iE "bam|bai|tabix|index" | tr '\n' ' ')
if [ -z "$STEPS" ]; then echo "  note  no step list in the paradoxa log, so unverified"
elif [ -z "$IDX$ISTEP" ]; then echo "  OK    no BAM/BAI/tabix step or file: the 2^29 limit does not touch scorpioides"
else echo "  note  CARP makes indexes ($ISTEP$IDX). If one fails on scorpioides it is a side step: the rule still accepts the run once the unified GFF and its statistics exist, and lists the failed step in carp.done"; fi

echo; echo "=== 5. dry run on all $(echo $SPECIES | wc -w) species (nothing is submitted)"
ALL=$(for sp in $SPECIES; do printf '%s ' "$REFS/${sp}_repeats/MANIFEST.txt"; done)
snakemake $ALL -n --nolock --rerun-triggers mtime --profile profiles/slurm > /tmp/carp_dry_$TS.log 2>&1
rc=$?
sed -n '/^Job stats:/,/^$/p' /tmp/carp_dry_$TS.log | sed 's/^/  /'
if [ $rc -ne 0 ]; then
  stop "dry run failed (exit $rc), log /tmp/carp_dry_$TS.log"
  grep -m1 -A4 -E "Exception|Error" /tmp/carp_dry_$TS.log | sed 's/^/    /'
else
  for sp in $SPECIES; do
    n=$(grep -cE "wildcards: species=${sp}[[:space:]]*$" /tmp/carp_dry_$TS.log)
    if [ -s "$REFS/${sp}_repeats/MANIFEST.txt" ]; then
      if [ "$n" -eq 0 ]; then echo "  OK    $sp: published, no jobs"
      else stop "$sp is published but would rerun ($n jobs): $(grep -B12 -E "wildcards: species=${sp}[[:space:]]*$" /tmp/carp_dry_$TS.log | grep -m1 'reason:' | cut -c1-150)"; fi
    elif [ "$n" -eq 2 ]; then echo "  OK    $sp: annotate + publish"
    else stop "$sp: dry run plans $n jobs, expected 2 (annotate + publish)"; fi
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
    snakemake --help 2>/dev/null | grep -q -- "--prioritize" && PRIO="--prioritize $REFS/Drosera_scorpioides_repeats/MANIFEST.txt";;
  esac
  setsid nohup bash -c "snakemake $NT --profile profiles/slurm --rerun-triggers mtime --keep-going --jobs 4 $PRIO; echo \"snakemake exit \$? at \$(date)\"" > "$LOG" 2>&1 < /dev/null &
  echo; echo "LAUNCHED:$NEW"
  echo "controller log: $RP/$LOG   checks: $REP"
  sleep 90
  "$RP/scripts/carp_status.sh"
fi
