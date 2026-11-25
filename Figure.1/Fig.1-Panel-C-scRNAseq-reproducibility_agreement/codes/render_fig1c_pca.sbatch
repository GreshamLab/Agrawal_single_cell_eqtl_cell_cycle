#!/bin/bash
#
#SBATCH --verbose
#SBATCH --job-name=render_fig1c_pca
#SBATCH --output=render_fig1c_pca.o%j
#SBATCH --error=render_fig1c_pca.e%j
#SBATCH --account=torch_pr_130_general
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb

# Re-knit the Fig 1C condition PCA Rmd. Submit from this directory:  sbatch render_fig1c_pca.sbatch
cd "${SLURM_SUBMIT_DIR}"
module load r/4.5.1
export R_LIBS_USER=/scratch/aa9679/R/x86_64-pc-linux-gnu/4.5
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8   # so the em dash in titles renders instead of <U+2014>

RMD="Fig1-PanelC-does-gene-experession-differ-bwn-condition.Rmd"
echo "Rendering ${RMD} on $(hostname) at $(date)"
Rscript -e "rmarkdown::render('${RMD}', knit_root_dir = getwd())"
rc=$?
echo "Finished ${RMD} at $(date) (exit ${rc})"
exit ${rc}
