#!/bin/bash
#SBATCH --job-name=render_Fig3
#SBATCH --output=/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle/Figure.3/render_all_Fig3_%j.log
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=96gb
#SBATCH --account=torch_pr_130_general
#SBATCH --mail-type=FAIL,END
# Re-renders Figure 3 (YPD.1) in order: classification (Step1) -> Panel A (Step2) -> Panel B (Step3).
# Stops at the first failure, since each step reads the tables of the one before.
set -e
F3=/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle/Figure.3
R="singularity exec --overlay /share/apps/images/r-4.6.1-ubuntu-26.04.sqf:ro /share/apps/images/ubuntu-26.04.sif /ext3/r/4.6.1/bin/Rscript"
A=$F3/Fig.3-Panel-A-C-D-cis-and-trans-eqtl-YPD-variant-type

echo "$(date '+%F %H:%M') START classification"
cd $A/codes && $R -e "rmarkdown::render('Fig3_panel_B_C_D_variant_eqtl_classification.Rmd')" > ../results/render_step1_YPD.1.log 2>&1
echo "$(date '+%F %H:%M') OK    classification"
echo "$(date '+%F %H:%M') START Panel A"
cd $A/codes && $R -e "rmarkdown::render('Fig.3-panel-A-cis_trans_eqtl_mapping_YPD.Rmd')" > ../results/render_step2_YPD.1.log 2>&1
echo "$(date '+%F %H:%M') OK    Panel A"
echo "$(date '+%F %H:%M') START Panel B"
cd $F3/Fig.3-Panel-B-raw_data_proof_YPD && $R -e "rmarkdown::render('Fig3-Panel-B-YPD.1_raw_eqtl_verification.Rmd')" > render.log 2>&1
echo "$(date '+%F %H:%M') OK    Panel B"
echo "$(date '+%F %H:%M') DONE"
