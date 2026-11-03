#!/bin/bash
#SBATCH --job-name=fig4_step1_theme
#SBATCH --account=torch_pr_130_general
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --output=logs/%x_%j.out

# Re-knit only Step1 (figures) of every Figure 4 dataset, e.g. after a plotting-theme change.
P=/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle
module purge
module load r/4.6.1
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
for ds in YPD.2 YPD_integrated Flu.1 Flu.2 Flu_integrated sc2glu.1_only sc2glu.2_only sc2glu_integrated scgal.1; do
  bash $P/shared_code/rerun_fig3_fig4/render_rmd.sh $P/Figure.4/Fig.4-${ds}/codes/Step1_Fig4_panel_B_C_D_variant_eqtl_classification_${ds}.Rmd
done
