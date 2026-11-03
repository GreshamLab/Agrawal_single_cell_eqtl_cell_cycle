#!/bin/sh
#
#SBATCH --verbose
#SBATCH --job-name=pileup_YPD.1
#SBATCH --output=pileup.o%j
#SBATCH --error=pileup.e%j
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=15gb
#SBATCH --array=0-214
#SBATCH --account=torch_pr_130_general

module purge
python3 1.parallel_pileup.py $SLURM_ARRAY_TASK_ID

