#!/bin/bash
#
#SBATCH --verbose
#SBATCH --job-name=render_fig5
#SBATCH --output=render_fig5.o%A_%a
#SBATCH --error=render_fig5.e%A_%a
#SBATCH --account=torch_pr_130_general
#SBATCH --array=0-1%1
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96gb
#SBATCH --mail-type=END
#SBATCH --mail-user=aa9679@nyu.edu

# Knit the Figure 5 Rmds (one Rmd per array task); part II reads the harmonic results, so tasks run one at a time (%1).
# Submit from the codes/ directory:  sbatch render_fig5_rmds.sbatch
cd "${SLURM_SUBMIT_DIR}"
module load r/4.5.1
# org.Sc.sgd.db and other user packages live here (same library the interactive session uses)
export R_LIBS_USER=/scratch/aa9679/R/x86_64-pc-linux-gnu/4.5
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8   # so the em dash in plot titles renders instead of <U+2014>

RMDS=(
  "Fig5.cell-cycle-harmonic.Rmd"
  "Fig5.-part-II-compare-with-ypd.2.Rmd"
)
RMD="${RMDS[$SLURM_ARRAY_TASK_ID]}"
echo "Rendering ${RMD} on $(hostname) at $(date)"

Rscript -e "rmarkdown::render('${RMD}', knit_root_dir = getwd())"
rc=$?
echo "Finished ${RMD} at $(date) (exit ${rc})"
exit ${rc}
