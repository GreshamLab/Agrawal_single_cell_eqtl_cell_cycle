#!/bin/bash
#SBATCH --verbose
#SBATCH --job-name=render_Fig3_BCD
#SBATCH --output=render_Fig3_BCD.o%j
#SBATCH --error=render_Fig3_BCD.e%j
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96gb
#SBATCH --time=4:00:00
#SBATCH --account=torch_pr_130_general
# Knit Fig3_panel_B_C_D_variant_eqtl_classification.Rmd (no titles/captions, "Gene CNVs" label)
set -e
cd /projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle/Figure.3/Fig.3-Panel-A-C-D-cis-and-trans-eqtl-YPD-variant-type/codes
R="singularity exec --overlay /share/apps/images/r-4.6.1-ubuntu-26.04.sqf:ro /share/apps/images/ubuntu-26.04.sif /ext3/r/4.6.1/bin/Rscript"
$R -e "rmarkdown::render('Fig3_panel_B_C_D_variant_eqtl_classification.Rmd')" > Fig3_panel_B_C_D_variant_eqtl_classification.render_2026-11-24_genecnv.log 2>&1
echo "DONE $(date)"
