#!/bin/sh
#
#SBATCH --job-name=fig4_best_rep
#SBATCH --output=../results/render_%x_%j.log
#SBATCH --error=../results/render_%x_%j.log
#SBATCH --account=torch_pr_130_general
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aa9679@nyu.edu

# Knit one best-replicate comparison Rmd (YPD.1, sc2glu.2, Flu.1, Flu.2). Submit from codes/:
#   sbatch --job-name=concordance render_best_replicates.sbatch Fig4.effect-size-comparison-best-replicates.Rmd
#   sbatch --job-name=partII      render_best_replicates.sbatch Fig4.PartII-compare-conditions-best-replicates.Rmd
cd "${SLURM_SUBMIT_DIR}"
module load r/4.6.1
Rscript -e "rmarkdown::render('$1')"
