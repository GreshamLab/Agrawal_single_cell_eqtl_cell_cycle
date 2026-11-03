#!/bin/sh
#
#SBATCH --verbose
#SBATCH --job-name=YPD1_Q4
#SBATCH --output=Q4_split.o%j
#SBATCH --error=Q4_split.e%j
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=5gb
#SBATCH --array=0-13

module purge
python3 split_bam.py $SLURM_ARRAY_TASK_ID
