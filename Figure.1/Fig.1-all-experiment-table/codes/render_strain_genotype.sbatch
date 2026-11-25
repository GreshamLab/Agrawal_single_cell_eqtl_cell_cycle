#!/bin/bash
#
#SBATCH --verbose
#SBATCH --job-name=render_strain_genotype
#SBATCH --output=render_strain_genotype.o%j
#SBATCH --error=render_strain_genotype.e%j
#SBATCH --account=torch_pr_130_general
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb

cd "${SLURM_SUBMIT_DIR}"
module load r/4.5.1
export R_LIBS_USER=/scratch/aa9679/R/x86_64-pc-linux-gnu/4.5
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
RMD="1a.Table_describing_the_strain_genotype.Rmd"
Rscript -e "rmarkdown::render('${RMD}', knit_root_dir = getwd(), output_dir = '/scratch/aa9679/TMPDIR/claude-3859618/-projects-rps-cgsb-gresham-Akriti-wildyeast1/41cd9467-c670-4998-b588-9042efc9db83/scratchpad')"
