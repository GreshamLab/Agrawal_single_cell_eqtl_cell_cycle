#!/bin/sh
#
#SBATCH --verbose
#SBATCH --job-name=demux_ypd2
#SBATCH --output=ypd2.o%j
#SBATCH --error=ypd2.e%j
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=35gb
#SBATCH --array=0-214
#SBATCH --account=torch_pr_130_general

module purge
python3 demuxlet.py $SLURM_ARRAY_TASK_ID

#find demuxlet_out_YPD.2 -name "demuxlet_output.tsv.best" -exec awk 'NR == 1 {print} $1 != "INT_ID" {print}' {} + > merged_demuxlet_output_YPD_2.tsv

#Another way:
#find demuxlet_out_YPD.2 -name "demuxlet_output.tsv.best" -exec awk 'FNR>1 || NR==1' {} + > merged_demuxlet_output_YPD_2_new.tsv