#!/bin/bash
set -euo pipefail
cd /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
eval "$(micromamba shell hook --shell bash)"
micromamba activate smk
echo "starting at $(date)"
snakemake --profile profiles/slurm --verbose --printshellcmds \
  results/Drosera_roseana/blobtoolkit/initial/hap1/.done
echo "finished at $(date)"
