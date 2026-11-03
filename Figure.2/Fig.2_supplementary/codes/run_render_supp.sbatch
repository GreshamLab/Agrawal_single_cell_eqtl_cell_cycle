#!/bin/bash
#SBATCH --job-name=render_fig2_supp
#SBATCH --account=torch_pr_130_general
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle/Figure.2/Fig.2_supplementary/codes/%x_%j.log

module purge
module load r/4.6.1
cd /projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle/Figure.2/Fig.2_supplementary/codes
Rscript -e "rmarkdown::render('fig.2.supplementary-replicates-per-condition-DEG.Rmd')"
