#!/bin/bash
# CARP repeat + rDNA annotation of D. paradoxa and D. binata (workflow/rules/repeats.smk)
# and its publication to refs/{species}_repeats/. One Slurm job per species.
# Run it in tmux on the login node; re-running resumes (CARP keeps its work):
#   tmux new-session -d -s carp 'bash scripts/launchers/launch_carp.sh; exec bash'
# Calls the smk environment's own snakemake: a shell that tmux starts has no
# micromamba shell hook, so 'micromamba activate' cannot be relied on here.
SMK_ENV=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk
export PATH="$SMK_ENV/bin:$PATH"
cd /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd || exit 1
REFS=/netscratch/dep_mercier/grp_marques/Aaryan/refs
mkdir -p logs
LOG=logs/launch_carp_$(date +%Y%m%d_%H%M).log
echo "started $(date) on $(hostname); $(command -v snakemake) $(snakemake --version)" | tee -a "$LOG"
# targets first: --allowed-rules takes every word after it as a rule name. It also
# keeps Snakemake from re-running anything upstream of the frozen assemblies.
snakemake "$REFS/Drosera_paradoxa_repeats/MANIFEST.txt" "$REFS/Drosera_binata_repeats/MANIFEST.txt" \
    --profile profiles/slurm --printshellcmds \
    --allowed-rules annotate_repeats_carp publish_repeats 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
echo "snakemake exit $rc at $(date)" | tee -a "$LOG"
