#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=10
#SBATCH --time=10:00:00
#SBATCH --mem=15GB
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%j.err
#SBATCH --job-name=YPD1_align
#SBATCH --account=torch_pr_130_general
  
module purge
module load cellranger/7.1.0
#make a new directory

#make a new directory
mkdir -p output_of_alignment_YPD1_with_clean_gtf
#move into that directory
cd output_of_alignment_YPD1_with_clean_gtf

#cellranger out to get the reference
cellranger count --id=run_count \
--transcriptome=/projects/rps/cgsb/gresham/Akriti/wildyeast1/creating_different_Ref_genomes_for_S288C/customized_ref_genome_with_kan_and_nat_annotations/S.cerevisiae_ref_with_clean_gtf \
--fastqs=/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/YPD.1/alignment/merged_files_YPD.1 \
--sample=H23N5DRX5_wildyeast_10X \
--localcores=10 \
--localmem=210

echo Run is finished
 
#took 8 hours to run and only 8GB with the above credential. If I set to this much memory and increase node it can be faster
