#!/bin/bash
#
#SBATCH --verbose
#SBATCH --job-name=render_figS456
#SBATCH --output=render_figS456.o%j
#SBATCH --error=render_figS456.e%j
#SBATCH --account=torch_pr_130_general
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb

# Re-knit the Fig S4-6 Rmd. Submit from this directory:  sbatch render_figS456.sbatch
cd "${SLURM_SUBMIT_DIR}"
module load r/4.5.1
export R_LIBS_USER=/scratch/aa9679/R/x86_64-pc-linux-gnu/4.5
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

RMD="Fig4.5.6.supplementary_strain_gene_expression_correlation.Rmd"
echo "Rendering ${RMD} on $(hostname) at $(date)"
Rscript -e "rmarkdown::render('${RMD}', knit_root_dir = getwd())"
rc=$?
echo "Finished ${RMD} at $(date) (exit ${rc})"
exit ${rc}
