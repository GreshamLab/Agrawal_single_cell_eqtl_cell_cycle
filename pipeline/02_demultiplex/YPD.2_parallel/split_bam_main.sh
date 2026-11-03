#!/bin/sh
#
#SBATCH --verbose
#SBATCH --job-name=split_bam
#SBATCH --output=split_bam.o%j
#SBATCH --error=split_bam.e%j
#SBATCH --time=1:30:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10gb
#SBATCH --array=0-214
#SBATCH --account=torch_pr_130_general


module purge
module load bioinformatics/20260224

python3 5.split_bam_new_code.py $SLURM_ARRAY_TASK_ID
