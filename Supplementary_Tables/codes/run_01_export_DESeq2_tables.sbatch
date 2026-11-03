#!/bin/bash
#SBATCH --job-name=export_DESeq2
#SBATCH --account=torch_pr_130_general
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=../logs/%x_%j.out

module purge
module load r/4.6.1
Rscript 01_export_DESeq2_tables.R
