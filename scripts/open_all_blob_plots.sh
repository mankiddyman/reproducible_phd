#!/bin/bash
# Open all blobtoolkit static plots (blob, snail, cumulative) in Firefox tabs.
# Run on a machine with the netscratch path mounted and firefox installed.
set -euo pipefail
cd /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd

# Collect plot paths grouped by species
plots=()
for sp_dir in results/Drosera_*; do
  sp=$(basename "$sp_dir")
  for hap_dir in "$sp_dir"/blobtoolkit/initial/*/output/blobtoolkit/plots; do
    [ -d "$hap_dir" ] || continue
    # Get hap1/hap2/p_utg from the path
    target=$(basename $(dirname $(dirname "$hap_dir")))
    for kind in blob snail cumulative; do
      plot="$hap_dir/$sp.$kind.png"
      [ -f "$plot" ] && plots+=("$plot")
    done
  done
done

if [ ${#plots[@]} -eq 0 ]; then
  echo "No plots found."
  exit 1
fi

echo "Opening ${#plots[@]} plots in Firefox..."
# Open all in one firefox call so they share a single window
firefox "${plots[@]}" &
