#!/bin/sh
#
#SBATCH --verbose
#SBATCH --job-name=pileup_YPD1_noSD
#SBATCH --output=split_bam_pileup.o%j
#SBATCH --error=split_bam_pileup.e%j
#SBATCH --time=14:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=4gb
#SBATCH --array=0-53


module purge
python3 parallel_pileup.py $SLURM_ARRAY_TASK_ID

