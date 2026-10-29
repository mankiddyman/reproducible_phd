#!/usr/bin/env bash
# carp_launch.sh -- CARP repeat + rDNA annotation (rule annotate_repeats_carp, then publish_repeats)
# for every species in config/annotation.csv not yet published to refs/. Checks first; submits only
# if every check passes:
#   0   this host resolves your user and Slurm accepts a submission from it (sbatch --test-only)
#   0b  every idle node the CARP jobs may use (a) can look up your user: hpc07 could not on 2026-10-28,
#       so Singularity died there in a second and the node attracted every job; and (b) leaves CARP a
#       memory budget of >= MIN_BUDGET GB. The nodes run strict overcommit (vm.overcommit_memory=2):
#       every fork is charged to the node's CommitLimit, part of which processes outside Slurm hold.
#       On 2026-10-28 DANTE's forks failed with Errno 12 on hpc09 (200 GB of room) and hpc08 (358 GB);
#       the per-species DANTE caps in config are sized to fit MIN_BUDGET
#   1   Snakemake submits through Slurm (profiles/slurm); nothing else holds the repo
#   2   the CARP image is the one the published species were annotated with (sha256 in their MANIFESTs)
#   3   each job's Slurm time limit is >= 1.5x its predicted runtime (D. paradoxa's rate, plus the extra
#       DANTE time the species' worker cap costs)
#   4   BAM/BAI/tabix in what CARP made for paradoxa (scorpioides has 2 chromosomes > 2^29 bp) -- informational
#   5   dry run restricted to the two CARP rules: published species get no jobs, every other species
#       exactly 2 (annotate + publish), and no other rule. The unrestricted dry run is reported alongside:
#       an edit to config/species.csv makes Snakemake want to rebuild every frozen assembly, and
#       --allowed-rules keeps CARP on the frozen FASTAs as they are
# Species whose CARP job is already in Slurm are not submitted again: the controller waits for those
# jobs to end, then runs the same two rules for them (publish if annotated, otherwise annotate again).
# One detached controller (setsid nohup), at most 4 species at a time, scorpioides first.
# Report: logs/carp_gate_<time>.txt. Status at any time: scripts/carp_status.sh
RP=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
ENV=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk
SIF=/netscratch/dep_mercier/grp_marques/Aaryan/methods/carp/carp_1.9.5.sif
REFS=/netscratch/dep_mercier/grp_marques/Aaryan/refs
CFG=config/config.yaml
CARP_RULES="annotate_repeats_carp publish_repeats"
MIN_BUDGET=250
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
Q=$(squeue -u "$ME" -h -o '%i %k' 2> /dev/null)
RUN_SP=""; RUN_IDS=""
for sp in $SPECIES; do
  id=$(echo "$Q" | awk -v s="$sp" '$2 ~ ("wildcards_" s "$") {print $1; exit}')
  [ -n "$id" ] && { RUN_SP="$RUN_SP $sp"; RUN_IDS="$RUN_IDS $id"; }
done

{
echo "CARP on every species: checks first, submit only if all pass   ($(date '+%F %T'), $(hostname))"
[ -n "$RUN_SP" ] && echo "  already in Slurm, followed up when their jobs end:$RUN_SP (job$RUN_IDS)"

echo; echo "=== 0. this host can run the controller"
if getent passwd "$ME" > /dev/null; then echo "  OK    $(hostname) resolves uid $ME"
else stop "$(hostname) cannot resolve your user (uid $ME). Run scripts/carp_launch.sh on hpc09"; fi
if sbatch --test-only -p normal -t 1 -c 1 --mem=100M --wrap true > /tmp/carp_sbtest_$TS.txt 2>&1; then
  echo "  OK    Slurm accepts a test submission from $(hostname)"
else stop "sbatch --test-only refused from $(hostname): $(head -1 /tmp/carp_sbtest_$TS.txt)"; fi

echo; echo "=== 0b. idle nodes CARP may use: user lookup, and CARP's memory budget there (>= $MIN_BUDGET GB)"
EXCL=$(python3 -c "import yaml; print((((yaml.safe_load(open('$CFG')) or {}).get('repeats') or {}).get('resources') or {}).get('default', {}).get('exclude', 'hpc15'))")
echo "  excluded for every CARP job: $EXCL   (? = no free slot within 20 s, not checked)"
PROBE='if getent passwd "$(id -u)" > /dev/null; then printf "NSS_OK "; else printf "NSS_FAIL "; fi; awk "/^CommitLimit:/{c=\$2} /^Committed_AS:/{a=\$2} /^MemAvailable:/{m=\$2} END{printf \"%d %d %d \", c, a, m}" /proc/meminfo; cat /proc/sys/vm/overcommit_memory'
BAD=""; LINE=""
for n in $(sinfo -h -p normal -t idle,mix -o '%N' | xargs -r -n1 scontrol show hostnames 2> /dev/null | sort -u); do
  case ",$EXCL," in *",$n,"*) LINE="$LINE $n:excluded"; continue;; esac
  r=$(timeout 45 srun -I20 -Q -w "$n" -p normal -t 1 -c 1 --mem=100M -J carp_nodecheck bash -c "$PROBE" 2> /dev/null)
  if [ -z "$r" ]; then LINE="$LINE $n:?"; continue; fi
  set -- $r
  b=$(awk -v c="$2" -v a="$3" -v m="$4" -v mode="$5" 'BEGIN { b = 810; room = (c - a) / 1048576; if (mode == 2 && 0.8 * room < b) b = 0.8 * room; if (0.9 * m / 1048576 < b) b = 0.9 * m / 1048576; printf "%d", b }')
  if [ "$1" != NSS_OK ]; then LINE="$LINE $n:NO-USER"; BAD="$BAD $n"
  elif [ "$b" -lt "$MIN_BUDGET" ]; then LINE="$LINE $n:${b}GB-LOW"; BAD="$BAD $n"
  else LINE="$LINE $n:${b}GB"; fi
done
echo "  $LINE" | fold -w 150 -s
if [ -n "$BAD" ]; then stop "idle nodes CARP may use but should not:$BAD -- add them to exclude under repeats.resources in $CFG"
else echo "  OK    every idle node CARP may use can look up your user and leaves CARP >= $MIN_BUDGET GB"; fi

echo; echo "=== 1. Snakemake submits to Slurm, nothing else holds the repo"
if grep -qE "^executor:[ \"']*slurm" profiles/slurm/config.yaml 2> /dev/null; then echo "  OK    profiles/slurm submits to Slurm"
else stop "profiles/slurm/config.yaml has no 'executor: slurm': jobs would run on $(hostname) itself"; fi
python3 -c "import snakemake_executor_plugin_slurm" 2> /dev/null && echo "  OK    slurm executor plugin, snakemake $(snakemake --version)" \
  || stop "snakemake-executor-plugin-slurm not importable from $ENV"
NL=$(ls .snakemake/locks 2> /dev/null | wc -l)
if [ "$NL" -gt 0 ]; then stop "$NL lock file(s) in .snakemake/locks: another Snakemake may be running in this repo ($(pgrep -u "$ME" -fc snakemake) snakemake processes on $(hostname); check the other nodes)"
else echo "  OK    no lock on the repo"; fi
if pgrep -u "$ME" -f "carp_controller_|snakemake.*_repeats/MANIFEST.txt" > /dev/null; then stop "a CARP controller is already running on $(hostname): see scripts/carp_status.sh"; fi

echo; echo "=== 2. same image as the published species were annotated with"
if [ -s "$SIF" ]; then
  IMG=$(sha256sum "$SIF" | cut -d' ' -f1)
  REC=$(grep -hoE "[0-9a-f]{64}" "$REFS"/*_repeats/MANIFEST.txt 2> /dev/null | sort -u)
  if [ -z "$REC" ]; then echo "  note  no sha256 in the MANIFESTs; image not compared ($IMG)"
  elif echo "$REC" | grep -qx "$IMG"; then echo "  OK    sha256 $IMG matches the MANIFESTs"
  else stop "image sha256 $IMG is not the one in the MANIFESTs ($(echo $REC))"; fi
else stop "image missing: $SIF"; fi

echo; echo "=== 3. time limits against the predicted runtime (DANTE cap included)"
python3 - "$CFG" "$REFS" << 'PY' || stop "a Slurm time limit is shorter than 1.5x the predicted runtime (table above)"
import csv, glob, os, re, sys, yaml
CFG, REFS = sys.argv[1], sys.argv[2]
res = ((yaml.safe_load(open(CFG)) or {}).get("repeats") or {}).get("resources") or {}
RATE = 9.2 / 4.30   # h per Gb: D. paradoxa on 2026-10-15 with 64 DANTE workers, the slower finished run
dh = 0.0
for b in glob.glob("results/Drosera_paradoxa/repeats/carp/benchmarks/*dante*.tsv"):
    rows = [l.rstrip("\n").split("\t") for l in open(b) if l.strip()]
    try: dh = max(dh, float(dict(zip(rows[0], rows[-1]))["s"]) / 3600)
    except Exception: pass
print(f"  paradoxa's longest DANTE step took {dh:.1f} h at 64 workers" if dh else
      "  no DANTE benchmark from paradoxa: the estimate assumes the whole run slows with the cap")
def hours(t):
    t, d = str(t), 0
    if "-" in t: d, t = t.split("-", 1); d = int(d)
    p = [int(x) for x in t.split(":")]
    if d: h = p[0] + (p[1] / 60 if len(p) > 1 else 0)
    elif len(p) == 3: h = p[0] + p[1] / 60
    else: h = p[0] / 60
    return 24 * d + h
bad = 0
print(f"  {'species':22s} {'state':9s} {'Gb':>5s} {'DANTE':>5s} {'~h':>6s} {'limit h':>8s}")
for row in csv.DictReader(open("config/annotation.csv")):
    sp = row["species_id"]
    r = dict(res.get("default") or {}); r.update(res.get(sp) or {})
    fai = f"results/{sp}/assembly_final/{row['frozen_stage']}/{sp}_chr.fa.fai"
    gb = sum(int(l.split("\t")[1]) for l in open(fai)) / 1e9 if os.path.exists(fai) else float("nan")
    done = os.path.isfile(f"{REFS}/{sp}_repeats/MANIFEST.txt")
    m = re.search(r"\bdante=(\d+)", str(r.get("carp_extra", "")))
    cap = int(m.group(1)) if m else 64
    base = RATE * gb
    est = base + (dh * gb / 4.30 * (64 / cap - 1) if dh else base * (64 / cap - 1))
    try: lim = hours(r["time"])
    except Exception: lim = float("nan")
    ok = done or (lim >= 1.5 * est)
    bad += not ok
    print(f"  {sp:22s} {'published' if done else 'to run':9s} {gb:5.2f} {cap:5d} {est:6.1f} {lim:8.1f}{'' if ok else '   <-- limit too short'}")
sys.exit(1 if bad else 0)
PY

echo; echo "=== 4. BAM/BAI/tabix in what CARP made for paradoxa (scorpioides: 2 chromosomes > 2^29 bp)"
STEPS=$(grep -ohE "^(local)?rule [A-Za-z0-9_]+" logs/annotate_repeats_carp/Drosera_paradoxa.log 2> /dev/null | awk '{print $2}' | sort -u)
IDX=$(find results/Drosera_paradoxa/repeats -type f \( -name "*.bai" -o -name "*.csi" -o -name "*.tbi" -o -name "*.bam" \) 2> /dev/null | head -5 | tr '\n' ' ')
ISTEP=$(echo "$STEPS" | grep -iE "bam|bai|tabix|index" | tr '\n' ' ')
if [ -z "$STEPS" ]; then echo "  note  no step list in the paradoxa log, so unverified"
elif [ -z "$IDX$ISTEP" ]; then echo "  OK    none of CARP's $(echo $STEPS | wc -w) steps or files is a BAM/BAI/tabix index"
else echo "  note  CARP makes indexes ($ISTEP$IDX); a failure there is a side step, which the rule tolerates"; fi

ALL=$(for sp in $SPECIES; do printf '%s ' "$REFS/${sp}_repeats/MANIFEST.txt"; done)
echo; echo "=== 5a. what an unrestricted run would also redo (reported only)"
snakemake $ALL -n --nolock --rerun-triggers mtime --profile profiles/slurm > logs/carp_dry_full_$TS.log 2>&1
jobs_tsv logs/carp_dry_full_$TS.log > /tmp/carp_jobs_full_$TS.tsv
echo "  $(wc -l < /tmp/carp_jobs_full_$TS.tsv) jobs, of which outside the CARP rules: $(awk -F'\t' '$1 != "annotate_repeats_carp" && $1 != "publish_repeats"' /tmp/carp_jobs_full_$TS.tsv | wc -l)"
awk -F'\t' '$1 != "annotate_repeats_carp" && $1 != "publish_repeats" && $3 !~ /^Input files updated by another job/ {print $1 ": " substr($3, 1, 150)}' \
  /tmp/carp_jobs_full_$TS.tsv | sort | uniq -c | sort -rn | head -6 | sed 's/^/   root cause: /'

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
    elif [ "$n" -eq 2 ]; then echo "  OK    $sp: annotate + publish, excluding ${ex:-nothing}"
    else stop "$sp: restricted dry run plans $n jobs, expected 2"; fi
  done
fi
} 2>&1 | tee "$REP"

NEW=""
for sp in $SPECIES; do
  [ -s "$REFS/${sp}_repeats/MANIFEST.txt" ] && continue
  case " $RUN_SP " in *" $sp "*) continue;; esac
  NEW="$NEW $sp"
done
if [ -s "$STOPF" ]; then
  echo; echo "NOT LAUNCHED: $(wc -l < "$STOPF") check(s) failed:"; sed 's/^/  - /' "$STOPF"
  echo "Paste back $REP"
elif [ -z "$NEW$RUN_SP" ]; then
  echo; echo "Nothing to launch: every species is published."
else
  LOG=logs/launch_carp_$TS.log
  CTRL=$RP/logs/carp_controller_$TS.sh
  NT=$(for sp in $NEW; do printf '%s ' "$REFS/${sp}_repeats/MANIFEST.txt"; done)
  RT=$(for sp in $RUN_SP; do printf '%s ' "$REFS/${sp}_repeats/MANIFEST.txt"; done)
  IDS=$(echo $RUN_IDS | tr ' ' ',')
  PRIO=""
  case " $NEW " in *" Drosera_scorpioides "*)
    snakemake --help 2> /dev/null | grep -q -- "--prioritize" && PRIO="--prioritize $REFS/Drosera_scorpioides_repeats/MANIFEST.txt";;
  esac
  cat > "$CTRL" << WEOF
#!/usr/bin/env bash
# CARP controller written by scripts/carp_launch.sh at $TS
export PATH=$ENV/bin:\$PATH
cd $RP || exit 1
if [ -n "$NT" ]; then
  snakemake $NT --allowed-rules $CARP_RULES --profile profiles/slurm --rerun-triggers mtime --keep-going --jobs 4 $PRIO
  echo "snakemake exit \$? at \$(date)"
fi
if [ -n "$RT" ]; then
  echo "waiting for the CARP jobs that were already running ($IDS)"
  while squeue -h -j "$IDS" 2> /dev/null | grep -q .; do sleep 600; done
  snakemake $RT --allowed-rules $CARP_RULES --profile profiles/slurm --rerun-triggers mtime --keep-going --jobs 4
  echo "follow-up snakemake exit \$? at \$(date)"
fi
WEOF
  setsid nohup bash "$CTRL" > "$LOG" 2>&1 < /dev/null &
  echo; echo "LAUNCHED:${NEW:- nothing new}"
  [ -n "$RUN_SP" ] && echo "followed up once their running job ends:$RUN_SP"
  echo "controller: $CTRL   log: $RP/$LOG   checks: $REP"
  sleep 90
  "$RP/scripts/carp_status.sh"
fi
