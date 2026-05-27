#!/bin/bash
set -euo pipefail

cd /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd

eval "$(micromamba shell hook --shell bash)"
micromamba activate smk

echo "snakemake: $(which snakemake)"
echo "snakemake version: $(snakemake --version)"
echo "cwd: $(pwd)"
echo "starting at $(date)"

# Only the 4 species that failed on hpc02
# binata, roseana, filiformis are already running as 788381/788380/788382
snakemake --profile profiles/slurm --verbose --printshellcmds \
  results/Drosera_aliciae/blobtoolkit/initial/hap1/.done \
  results/Drosera_aliciae/blobtoolkit/initial/hap2/.done \
  results/Drosera_paradoxa/blobtoolkit/initial/hap1/.done \
  results/Drosera_paradoxa/blobtoolkit/initial/hap2/.done \
  results/Drosera_scorpioides/blobtoolkit/initial/hap1/.done \
  results/Drosera_scorpioides/blobtoolkit/initial/hap2/.done \
  results/Drosera_tokaiensis/blobtoolkit/initial/hap1/.done \
  results/Drosera_tokaiensis/blobtoolkit/initial/hap2/.done

echo "finished at $(date)"
