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

# Target: the all rule in workflow/Snakefile, which uses ploidy-aware target
# selection (decontamination_targets per species). Was previously a hardcoded
# list of hap1/hap2 across all species, which spawned meaningless jobs for
# polyploid species (scorpioides, aliciae, tokaiensis).
snakemake --profile profiles/slurm --verbose --printshellcmds all

echo "finished at $(date)"
