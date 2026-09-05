#!/bin/bash
#SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --time=00:30:00
#SBATCH --mem=30GB
#SBATCH --job-name=inferelator_test
#SBATCH --account=torch_pr_130_general
#SBATCH --output=slurm_program_times.o
#SBATCH --error=slurm_program_times.e
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=aa9679@nyu.edu

module purge

singularity exec --fakeroot \
    --overlay /scratch/$USER/my_env/overlay-15GB-500K.ext3:rw \
    /share/apps/images/cuda12.3.2-cudnn9.0.0-ubuntu-22.04.4.sif \
    /bin/bash -c "source /ext3/env.sh && conda activate inferelator_env && python 1.1_cc_time_all_genes_with_all_spellman_genes.py"



# singularity exec --fakeroot \
#     --overlay /scratch/$USER/my_env/overlay-15GB-500K.ext3:rw \
#     /share/apps/images/cuda12.3.2-cudnn9.0.0-ubuntu-22.04.4.sif \
#     /bin/bash -c "source /ext3/env.sh && conda activate inferelator_env && python export_all_gene_programs.py"
