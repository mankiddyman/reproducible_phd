#!/bin/bash
#SBATCH --job-name=hic_lib_joins
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=logs/hic_library_joins_%j.out
# chr1/chr2 joins and chr5 x chr6, per Hi-C library, from the re-mapped pairs of the weekend job
cd /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk/bin/python workflow/scripts/hic_library_joins.py \
    results/Drosera_paradoxa/qc/hic_remap/work/pairs_all.tsv results/Drosera_paradoxa/qc/hic_remap/library_joins.txt
