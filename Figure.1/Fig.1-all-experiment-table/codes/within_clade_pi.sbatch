#!/bin/bash
#
#SBATCH --verbose
#SBATCH --job-name=within_clade_pi
#SBATCH --output=within_clade_pi.o%j
#SBATCH --error=within_clade_pi.e%j
#SBATCH --account=torch_pr_130_general
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48gb

# Pairwise nucleotide-diversity matrices (no MAF filter) from the 96-strain all-sites VCF.
# Output is read by the within-clade pi chunk of 1a.Table_describing_the_strain_genotype.Rmd.
cd "${SLURM_SUBMIT_DIR}"
module load anaconda3/2025.06
python3 within_clade_pi.py \
  /projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/vcf_files_sent_by_Anne/96strains.vcf \
  ../results/within_clade_pi_matrices 16
