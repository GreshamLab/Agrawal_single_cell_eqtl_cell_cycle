#!/bin/bash
#
#SBATCH --verbose
#SBATCH --job-name=render_fig1_supp
#SBATCH --output=render_fig1_supp.o%j
#SBATCH --error=render_fig1_supp.e%j
#SBATCH --account=torch_pr_130_general
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb

# Re-knit the Fig 1 supplementary Rmds. Submit from this directory:  sbatch render_Fig1_supplementary.sbatch [file.Rmd ...]
# With no arguments both Rmds are rendered.
cd "${SLURM_SUBMIT_DIR}"
module load r/4.5.1
export R_LIBS_USER=/scratch/aa9679/R/x86_64-pc-linux-gnu/4.5
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8   # so the em dash in titles renders instead of <U+2014>

status=0
RMDS=("$@")
[ ${#RMDS[@]} -eq 0 ] && RMDS=(Fig1.2.3.supplementary_demuxlet.Rmd
                              Fig4.5.6.supplementary_strain_gene_expression_correlation.Rmd)
for RMD in "${RMDS[@]}"; do
  echo "Rendering ${RMD} on $(hostname) at $(date)"
  Rscript -e "rmarkdown::render('${RMD}', knit_root_dir = getwd())"
  rc=$?
  echo "Finished ${RMD} at $(date) (exit ${rc})"
  [ ${rc} -ne 0 ] && status=${rc}
done
exit ${status}
