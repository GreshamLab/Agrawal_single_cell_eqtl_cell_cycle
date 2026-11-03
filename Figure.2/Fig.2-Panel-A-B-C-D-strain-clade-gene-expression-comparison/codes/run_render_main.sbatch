#!/bin/bash
#SBATCH --job-name=render_fig2_main
#SBATCH --account=torch_pr_130_general
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle/Figure.2/Fig.2-Panel-A-B-C-D-strain-clade-gene-expression-comparison/codes/%x_%j.log

module purge
module load r/4.6.1
cd /projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle/Figure.2/Fig.2-Panel-A-B-C-D-strain-clade-gene-expression-comparison/codes
Rscript -e "rmarkdown::render('Fig.2-panel-A-B-C-D-differential-gene-expression.Rmd')"
