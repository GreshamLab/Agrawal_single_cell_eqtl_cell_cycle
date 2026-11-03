#!/bin/bash
#
#SBATCH --job-name=pb_nofilt
#SBATCH --output=../results/logs_pseudobulk_no_cell_filter/%x_%A_%a.out
#SBATCH --error=../results/logs_pseudobulk_no_cell_filter/%x_%A_%a.err
#SBATCH --account=torch_pr_130_general
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64gb
#SBATCH --array=0-11

# Tasks 0-5: min_cells = 0 (no cell filter) -> new file next to the dataset's pseudobulk files
# Tasks 6-11: min_cells = 5 (validation) -> must reproduce the existing filtered matrix
# Submit from codes/:  mkdir -p ../results/logs_pseudobulk_no_cell_filter && sbatch 0_pseudobulk_no_cell_filter.sbatch

BASE=/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all
DATASETS=(YPD.1 YPD.2 Flu.1 Flu.2 scgal.1 sc2glu_integrated)
OUTS=(
  $BASE/YPD.1/fastlmm/data/pheno_data/pseudobulk_YPD.1_all_strains_no_cell_filter.txt
  $BASE/YPD.2/fastlmm/data/pheno_data/pseudobulk_YPD.2_all_strains_no_cell_filter.txt
  $BASE/Flu.1/fastlmm/data/pheno_data/pseudobulk_flu.1_all_strains_no_cell_filter.txt
  $BASE/Flu.2/fastlmm/data/pheno_data/pseudobulk_flu.2_all_strains_no_cell_filter.txt
  $BASE/scgal.1/fastlmm/data/pheno_data/pseudobulk_scgal.1_all_strains_no_cell_filter.txt
  $BASE/sc2glu.2/fastlmm/data/pheno_data/pseudobulked_sc2glu_integrated_all_strains_no_cell_filter.csv
)
VALIDATION_DIR=../results/validation_pseudobulk_min_cells_5

i=$((SLURM_ARRAY_TASK_ID % 6))
if [ "$SLURM_ARRAY_TASK_ID" -lt 6 ]; then
  MIN=0; OUT=${OUTS[$i]}
else
  MIN=5; mkdir -p $VALIDATION_DIR; OUT=$(realpath $VALIDATION_DIR)/$(basename ${OUTS[$i]})
fi

module purge
module load r/4.6.1
cd "${SLURM_SUBMIT_DIR}"
Rscript 0_pseudobulk_no_cell_filter.R ${DATASETS[$i]} $MIN $OUT
