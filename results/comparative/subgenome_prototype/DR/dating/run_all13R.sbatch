#!/bin/bash
#SBATCH --job-name=all13R
#SBATCH --array=1-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype/DR/dating/all13R_run%a/slurm.log
D=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/subgenome_prototype/DR/dating/all13R_run${SLURM_ARRAY_TASK_ID}
cd "$D" || exit 1
echo "host $(hostname)  start $(date)  dir $D"
/usr/bin/mcmctree mcmctree.ctl > run.log 2>&1
echo "exit $?  end $(date)"
ls -la mcmc.txt 2>/dev/null || echo "NO mcmc.txt"
