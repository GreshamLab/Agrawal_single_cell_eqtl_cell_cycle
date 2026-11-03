#!/bin/bash
#SBATCH --job-name=build_supp_xlsx
#SBATCH --account=torch_pr_130_general
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=../logs/%x_%j.out

module purge
module load r/4.6.1
Rscript 02_build_supplementary_workbook.R
