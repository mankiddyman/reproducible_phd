#!/bin/bash
set -euo pipefail

cd /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd

# Activate micromamba env so snakemake is on PATH
eval "$(micromamba shell hook --shell bash)"
micromamba activate smk

# Verify env is good
echo "snakemake: $(which snakemake)"
echo "snakemake version: $(snakemake --version)"
echo "cwd: $(pwd)"
echo "starting at $(date)"

# Targets: 14 sentinels (7 species × 2 haplotypes)
snakemake --profile profiles/slurm --verbose --printshellcmds \
  results/Drosera_binata/blobtoolkit/initial/hap1/.done \
  results/Drosera_binata/blobtoolkit/initial/hap2/.done \
  results/Drosera_paradoxa/blobtoolkit/initial/hap1/.done \
  results/Drosera_paradoxa/blobtoolkit/initial/hap2/.done \
  results/Drosera_roseana/blobtoolkit/initial/hap1/.done \
  results/Drosera_roseana/blobtoolkit/initial/hap2/.done \
  results/Drosera_scorpioides/blobtoolkit/initial/hap1/.done \
  results/Drosera_scorpioides/blobtoolkit/initial/hap2/.done \
  results/Drosera_aliciae/blobtoolkit/initial/hap1/.done \
  results/Drosera_aliciae/blobtoolkit/initial/hap2/.done \
  results/Drosera_tokaiensis/blobtoolkit/initial/hap1/.done \
  results/Drosera_tokaiensis/blobtoolkit/initial/hap2/.done \
  results/Drosera_filiformis/blobtoolkit/initial/hap1/.done \
  results/Drosera_filiformis/blobtoolkit/initial/hap2/.done

echo "finished at $(date)"
