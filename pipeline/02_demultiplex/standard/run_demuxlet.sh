#!/bin/sh
#
#SBATCH --verbose
#SBATCH --job-name=YPD1demuxlet1sd
#SBATCH --output=1sd.o%j
#SBATCH --error=1sd.e%j
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=12gb
#SBATCH --array=0-38

module purge
python3 demuxlet.py $SLURM_ARRAY_TASK_ID


#find demuxlet_YPD1_out_noSD_default -name "demuxlet_output.tsv.best" -exec awk 'NR == 1 {print} $1 != "INT_ID" {print}' {} + > merged_demuxlet_output_YPD_2.tsv
#Used 9gb for 7hours
