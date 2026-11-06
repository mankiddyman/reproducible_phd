#!/usr/bin/env bash
# expansion_status.sh -- where the 2026-11 expansion runs stand (stage A on gpu01: annotation + FANTASIA;
# stage B on hpc09/Slurm: OrthoFinder on 19 species). Reads logs/expansion/ only; takes a second.
#   bash scripts/launchers/expansion_status.sh
REPO=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
RUN=$REPO/logs/expansion
cd "$REPO" || exit 1
NEW="Luzula_sylvatica Juncus_effusus Cuscuta_epithymum Cuscuta_campestris Chionographis_japonica Chamaelirium_luteum Morus_notabilis Morus_alba"

echo "== markers"
for m in A.preflight A1.ok A1.failed A2.ok A2.failed B.blocked B.ok B.failed; do
  [ -e "$RUN/$m" ] && { echo "-- $m ($(date -r "$RUN/$m" '+%F %T'))"; sed 's/^/   /' "$RUN/$m" | head -n 8; }
done
echo "== progress (last 'steps done' line of each Snakemake run)"
for l in A1 A2 B; do
  f=$RUN/$l.snakemake.log
  [ -s "$f" ] && echo "   $l: $(grep -E '[0-9]+ of [0-9]+ steps' "$f" | tail -n 1)   errors: $(grep -c '^Error' "$f")   (updated $(date -r "$f" '+%F %T'))"
done
echo "== per species: final proteome / FANTASIA"
for s in $NEW; do
  p=results/$s/annotation/final/$s.final.proteins.fa
  g=results/$s/annotation/function/fantasia/$s.results.csv
  printf "   %-24s %s   %s\n" "$s" "$([ -s "$p" ] && echo "proteins $(grep -c '^>' "$p")" || echo 'proteins -')" "$([ -s "$g" ] && echo 'FANTASIA done' || echo 'FANTASIA -')"
done
echo "== Slurm (OrthoFinder)"
squeue -u "$USER" -h -o '   %i %k %T %M %R' 2> /dev/null | grep -i orthofinder || echo "   no OrthoFinder job in the queue"   # %k: the Slurm comment, rule_<name>_wildcards_...
for l in trees run; do   # the two OrthoFinder jobs: orthofinder_trees, then orthofinder
  [ -s logs/orthofinder/$l.log ] && echo "   OrthoFinder $l.log: $(grep -v '^\s*$' logs/orthofinder/$l.log | tail -n 1 | cut -c1-150)"
done
