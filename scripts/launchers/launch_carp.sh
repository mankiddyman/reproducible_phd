#!/bin/bash
# CARP repeat + rDNA annotation of D. paradoxa and D. binata (workflow/rules/repeats.smk)
# and its publication to refs/{species}_repeats/. One Slurm job per species.
# Run it in tmux on the login node; re-running resumes (CARP keeps its work):
#   tmux new-session -d -s carp 'bash scripts/launchers/launch_carp.sh; exec bash'
cd /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd || exit 1
eval "$(micromamba shell hook --shell bash)"
micromamba activate smk
REFS=/netscratch/dep_mercier/grp_marques/Aaryan/refs
mkdir -p logs
LOG=logs/launch_carp_$(date +%Y%m%d_%H%M).log
echo "started $(date) on $(hostname); snakemake $(snakemake --version)" | tee -a "$LOG"
# targets first: --allowed-rules takes every word after it as a rule name. It also
# keeps Snakemake from re-running anything upstream of the frozen assemblies.
snakemake "$REFS/Drosera_paradoxa_repeats/MANIFEST.txt" "$REFS/Drosera_binata_repeats/MANIFEST.txt" \
    --profile profiles/slurm --printshellcmds \
    --allowed-rules annotate_repeats_carp publish_repeats 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
echo "snakemake exit $rc at $(date)" | tee -a "$LOG"
