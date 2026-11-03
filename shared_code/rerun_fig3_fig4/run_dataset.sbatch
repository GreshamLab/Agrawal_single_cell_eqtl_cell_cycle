#!/bin/bash
#SBATCH --account=torch_pr_130_general
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96gb
#SBATCH --output=logs/%x_%j.out

# Per-dataset eQTL pipeline: Step1 (clumping, local/distant, summaries) -> Step2 (panel A
# maps) -> Step3 (raw single-cell verification). Usage: sbatch -J <dataset> run_dataset.sbatch <dataset>
P=/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle
ds="$1"
module purge
module load r/4.6.1
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

if [ "$ds" = "YPD.1" ]; then
  C=$P/Figure.3/Fig.3-Panel-A-C-D-cis-and-trans-eqtl-YPD-variant-type/codes
  steps=("$C/Fig3_panel_B_C_D_variant_eqtl_classification.Rmd"
         "$C/Fig.3-panel-A-cis_trans_eqtl_mapping_YPD.Rmd"
         "$P/Figure.3/Fig.3-Panel-B-raw_data_proof_YPD/Fig3-Panel-B-YPD.1_raw_eqtl_verification.Rmd")
else
  C=$P/Figure.4/Fig.4-${ds}/codes
  steps=("$C/Step1_Fig4_panel_B_C_D_variant_eqtl_classification_${ds}.Rmd"
         "$C/Step2_Fig4-panel-A-cis_trans_eqtl_mapping_${ds}.Rmd"
         "$C/Step3_Fig4-Panel-B-${ds}_raw_eqtl_verification.Rmd")
fi

status=0
for s in "${steps[@]}"; do
  bash $P/shared_code/rerun_fig3_fig4/render_rmd.sh "$s" || { status=1; [[ "$s" == *Step1* || "$s" == *classification.Rmd ]] && break; }
done
echo "$ds: $([ $status -eq 0 ] && echo ALL DONE || echo 'SOME STEP FAILED')"
exit $status
