#!/bin/bash
#
#SBATCH --verbose
#SBATCH --job-name=step3_sc2glu2
#SBATCH --output=step3_sc2glu2.o%j
#SBATCH --error=step3_sc2glu2.e%j
#SBATCH --account=torch_pr_130_general
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --mail-type=END
#SBATCH --mail-user=aa9679@nyu.edu

# Render the sc2glu.2-only bulk vs single-cell effect-size comparison.
# Submit from this directory:  sbatch render_step3_sc2glu2_only.sbatch
cd "${SLURM_SUBMIT_DIR}"
module load r/4.5.1
export R_LIBS_USER=/scratch/aa9679/R/x86_64-pc-linux-gnu/4.5
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8   # so the em dash in plot titles renders instead of <U+2014>

RMD="Step3.effect-size-correlation-bulk-vs-sc2glu.2-only.Rmd"
echo "Rendering ${RMD} on $(hostname) at $(date)"
Rscript -e "rmarkdown::render('${RMD}', knit_root_dir = getwd())"
rc=$?
echo "Finished ${RMD} at $(date) (exit ${rc})"
exit ${rc}
