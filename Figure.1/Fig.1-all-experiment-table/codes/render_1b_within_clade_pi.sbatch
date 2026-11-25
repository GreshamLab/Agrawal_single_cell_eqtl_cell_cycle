#!/bin/bash
#
#SBATCH --verbose
#SBATCH --job-name=render_1b_pi
#SBATCH --output=render_1b_pi.o%j
#SBATCH --error=render_1b_pi.e%j
#SBATCH --account=torch_pr_130_general
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128gb

# Knit 1b.Within_clade_SNP_diversity.Rmd (within-clade SNP pi from the unfiltered 96-strain VCF).
cd "${SLURM_SUBMIT_DIR}"
module load r/4.5.1
export R_LIBS_USER=/scratch/aa9679/R/x86_64-pc-linux-gnu/4.5
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
Rscript -e "rmarkdown::render('1b.Within_clade_SNP_diversity.Rmd', knit_root_dir = getwd())"
