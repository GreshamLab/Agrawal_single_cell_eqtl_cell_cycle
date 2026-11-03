#!/bin/sh
#
#SBATCH --job-name=fig4_cond
#SBATCH --output=logs_pipeline/fig4_%x_%j.out
#SBATCH --error=logs_pipeline/fig4_%x_%j.err
#SBATCH --account=torch_pr_130_general
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aa9679@nyu.edu

# One condition per job (Step1 -> PLINK clumping -> Step1/2/3 knits).
# Submit from Figure.4/:
#   mkdir -p logs_pipeline
#   for c in Flu.1 Flu.2 Flu_integrated sc2glu_integrated sc2glu.1_only sc2glu.2_only scgal.1 YPD.2 YPD_integrated; do
#     sbatch --job-name=$c run_condition_eqtl_pipeline.sbatch $c
#   done
cd "${SLURM_SUBMIT_DIR}"
bash run_condition_eqtl_pipeline.sh "$1"
