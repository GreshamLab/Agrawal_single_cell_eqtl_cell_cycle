#!/usr/bin/env python3

#Use --double-id \ because strains in the PLINK matrix share the same FID and IID (e.g., SACE_GAL   SACE_GAL) for all other than SNPs
#Don't use --double-id \ because the strain name is split at the underscore in the plink matrix (e.g., SACE    GAL) for SNPs - don't
import os

def create_batch_scripts_for_all_variants():
    """Create GWAS batch processing scripts for all variant types"""
    
    # Configuration - variant types with genotype paths and output directories
    variant_types = {
        'SNPs': {
            'genotype': '/projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/pan_transcriptome_paper/SNPs/SNPs_no_change_name_and_no_strain_filtering_only_MAF15_filtered/full969Matrix.SNPs.Biallelic.Var.DP10.99percNonMissing.MAF15.NoSubTelo20kb.plink',
            'outdir': '/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/YPD.1/fastlmm/results/results_SNPs_968_with_double_id'
        },
        'CNVs': {
            'genotype': '/projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/pan_transcriptome_paper/CNVs/CNVs_no_change_name_and_no_strain_filtering_only_MAF15_filtered/CNV_969strains_gain_maf15_N400_20221121.plink',
            'outdir': '/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/YPD.1/fastlmm/results/results_CNVs_with_double_id'
        },
        'Indels': {
            'genotype': '/projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/vcf_files_sent_by_Anne/SNPs_indels/modified_snp_indel_plink_for_96strains/Indels_only_96/indels.96samples',
            'outdir': '/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/YPD.1/fastlmm/results/results_Indels_with_double_id'
        },
        'SVs': {
            'genotype': '/projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/vcf_files_sent_by_Anne/SVs_from_SV.1087/SV.1087samples_updated_fam_files/SV.1087samples.SVs',
            'outdir': '/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/YPD.1/fastlmm/results/results_SVs_with_double_id'
        }
    }
    
    # Phenotype batch directory
    batch_dir = "/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/YPD.1/fastlmm/data/pheno_batches"
    script_output_dir = "/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/YPD.1/fastlmm/codes/batch_scripts_with_double_id"
    
    # Create output directory
    os.makedirs(script_output_dir, exist_ok=True)
    
    
    print("CREATING GWAS BATCH SCRIPTS FOR ALL VARIANT TYPES")
    
    
    # Count batch files
    batch_files = [f for f in os.listdir(batch_dir) 
                   if f.startswith("phenotype_batch_") and f.endswith(".txt")
                   and "summary" not in f]
    num_batches = len(batch_files)
    
    print(f"\nBatch directory: {batch_dir}")
    print(f"Found {num_batches} batch files")
    print(f"Found {len(variant_types)} variant types: {', '.join(variant_types.keys())}\n")
    
    if num_batches == 0:
        print("ERROR: No batch files found!")
        return
    
    job_info = []
    
    # Iterate over variant types
    for variant_type, variant_config in variant_types.items():
        print(f"\n{'=' * 80}")
        print(f"Processing variant type: {variant_type}")
        print(f"{'=' * 80}")
        print(f"Genotype: {variant_config['genotype']}")
        print(f"Output directory: {variant_config['outdir']}\n")
        
        genotype_path = variant_config['genotype']
        output_dir = variant_config['outdir']
        
        # Create run-task.bash script
        run_task_script = os.path.join(script_output_dir, f"run_task_{variant_type}.bash")
        
        run_task_content = f"""#!/bin/bash
export PST_NUM_THREADS=3

# Set variables
GENOTYPE="{genotype_path}"
OUTDIR="{output_dir}"
NPERM=100
THREADS=3

# Batch file directory
BATCH_DIR="{batch_dir}"
BATCH_FILE="${{BATCH_DIR}}/phenotype_batch_$(printf '%02d' $SLURM_ARRAY_TASK_ID).txt"

# Create output directory
mkdir -p $OUTDIR

echo "=== GWAS BATCH JOB $SLURM_ARRAY_TASK_ID FOR {variant_type} ===";
echo "Variant type: {variant_type}"
echo "Batch file: $BATCH_FILE"

# Check if batch file exists
if [ ! -f "$BATCH_FILE" ]; then
    echo "ERROR: Batch file $BATCH_FILE not found!"
    echo "Available batch files:"
    ls -la $BATCH_DIR/
    exit 1
fi

# Count files in this batch
BATCH_SIZE=$(wc -l < "$BATCH_FILE")
echo "Processing $BATCH_SIZE phenotype files in this batch"
echo "Expected duration: ~$((BATCH_SIZE * 87 / 60)) minutes"

# Read each phenotype file from the batch and process it
COUNTER=1
while IFS= read -r PHENOTYPE; do
    # Remove any trailing whitespace/newlines
    PHENOTYPE=$(echo "$PHENOTYPE" | tr -d '\\r\\n')
    
    # Skip empty lines
    if [ -z "$PHENOTYPE" ]; then
        continue
    fi
    
    GENE_NAME=$(basename "$PHENOTYPE" .norm.phen)
    
    echo
    echo "--- Processing file $COUNTER/$BATCH_SIZE: $GENE_NAME ---"
    echo "Phenotype file: $PHENOTYPE"
    
    # Check if file exists
    if [ ! -f "$PHENOTYPE" ]; then
        echo "ERROR: Phenotype file $PHENOTYPE does not exist!"
        COUNTER=$((COUNTER + 1))
        continue
    fi
    
    # Run GWAS for this phenotype
    echo "Running GWAS for $GENE_NAME..."
    START_TIME=$(date +%s)
    
    python /scratch/aa9679/victor-material/Caudal2024_repro/runFaSTLMM.v4.0.py \\
        -g $GENOTYPE \\
        -p $PHENOTYPE \\
        -o $OUTDIR \\
        --double-id \\
        --nb-permutations $NPERM \\
        --uncompressed \\
        --threads $THREADS
    
    EXIT_CODE=$?
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    
    if [ $EXIT_CODE -eq 0 ]; then
        echo "SUCCESS: $GENE_NAME completed in ${{DURATION}}s"
    else
        echo "FAILED: $GENE_NAME failed with exit code $EXIT_CODE after ${{DURATION}}s"
    fi
    
    COUNTER=$((COUNTER + 1))
    
done < "$BATCH_FILE"

echo
echo "=== BATCH $SLURM_ARRAY_TASK_ID COMPLETED FOR {variant_type} ==="
echo "Finished at: $(date)"
echo "Processed $((COUNTER - 1)) phenotype files"

# Create a completion marker
touch "$OUTDIR/batch_${{SLURM_ARRAY_TASK_ID}}_completed.marker"
"""
        
        with open(run_task_script, 'w') as f:
            f.write(run_task_content)
        
        os.chmod(run_task_script, 0o755)
        print(f"Created: {os.path.basename(run_task_script)}")
        
        # Create SLURM submission script
        slurm_script = os.path.join(script_output_dir, f"submit_gwas_{variant_type}.sh")
        
        slurm_content = f"""#!/bin/bash
#
#SBATCH --verbose
#SBATCH --job-name=GWAS_{variant_type}
#SBATCH --output=gwas_{variant_type}_%A_%a.o
#SBATCH --error=gwas_{variant_type}_%A_%a.e
#SBATCH --time=20:00:00       
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=8gb
#SBATCH --account=torch_pr_130_general
#SBATCH --mail-type=FAIL,END
#SBATCH --array=1-{num_batches}

# Load required modules
module purge

# Run the task script
/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/singularity/run-fastlmm-0.6.12.bash bash {run_task_script}
"""
        
        with open(slurm_script, 'w') as f:
            f.write(slurm_content)
        
        os.chmod(slurm_script, 0o755)
        print(f"Created: {os.path.basename(slurm_script)}")
        
        job_info.append({
            'variant_type': variant_type,
            'batches': num_batches,
            'submit_script': slurm_script
        })
    
    # Create master submission script
    master_script = os.path.join(script_output_dir, "submit_all_variants.sh")
    
    master_content = """#!/bin/bash
# Master script to submit GWAS jobs for all variant types

echo "======================================================================"
echo "SUBMITTING GWAS JOBS FOR ALL VARIANT TYPES"
echo "======================================================================"
echo

"""
    
    for info in sorted(job_info, key=lambda x: x['variant_type']):
        master_content += f"""echo "Submitting {info['variant_type']} ({info['batches']} batches)..."
sbatch {info['submit_script']}
sleep 1

"""
    
    master_content += """echo
echo "======================================================================"
echo "ALL JOBS SUBMITTED"
echo "======================================================================"
echo "Monitor jobs with: squeue -u $USER"
echo "Check results in corresponding results_[VARIANT_TYPE]/"
"""
    
    with open(master_script, 'w') as f:
        f.write(master_content)
    
    os.chmod(master_script, 0o755)
    
    # Create README
    readme_file = os.path.join(script_output_dir, "README.txt")
    readme_content = f"""GWAS Batch Processing Scripts for Variant Types
===========================================================================

Generated scripts for {len(variant_types)} variant types
Total batches: {num_batches}
Total jobs: {len(job_info) * num_batches}

VARIANT TYPES:
"""
    
    for info in sorted(job_info, key=lambda x: x['variant_type']):
        variant_type = info['variant_type']
        config = variant_types[variant_type]
        readme_content += f"\n{variant_type}:"
        readme_content += f"\n  Genotype: {config['genotype']}"
        readme_content += f"\n  Output: {config['outdir']}"
        readme_content += f"\n  Batches: {info['batches']}"
    
    readme_content += f"""

DIRECTORY STRUCTURE:
Each variant type has:
  - submit_gwas_[VARIANT].sh (SLURM submission script)
  - run_task_[VARIANT].bash (processing script)

USAGE:

1. Submit ALL variant types:
   bash submit_all_variants.sh

2. Submit individual variant type:
   sbatch submit_gwas_SNPs.sh
   sbatch submit_gwas_Indels.sh
   sbatch submit_gwas_CNV.sh
   sbatch submit_gwas_SVs.sh

3. Monitor jobs:
   squeue -u $USER

4. Check results:
   ls -la results_[VARIANT_TYPE]/

OUTPUT LOCATIONS:


PARAMETERS:
- Permutations: 100
- Threads: 1
- Memory: 15GB per job
- Time limit: 6 hours
- Array jobs: {num_batches} batches per variant type
"""
    
    with open(readme_file, 'w') as f:
        f.write(readme_content)
    
    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Created scripts for {len(variant_types)} variant types\n")
    
    for info in sorted(job_info, key=lambda x: x['variant_type']):
        print(f"{info['variant_type']:10s}: {info['batches']:3d} batches")
    
    print(f"\n{'=' * 80}")
    print(f"All scripts saved to: {script_output_dir}")
    print(f"{'=' * 80}\n")
    
    print("SUBMISSION OPTIONS:\n")
    print("1. Submit ALL variant types at once:")
    print(f"   bash submit_all_variants.sh\n")
    
    print("2. Submit individual variant type:")
    for info in sorted(job_info, key=lambda x: x['variant_type']):
        print(f"   sbatch submit_gwas_{info['variant_type']}.sh")
    print()

if __name__ == "__main__":
    create_batch_scripts_for_all_variants()
