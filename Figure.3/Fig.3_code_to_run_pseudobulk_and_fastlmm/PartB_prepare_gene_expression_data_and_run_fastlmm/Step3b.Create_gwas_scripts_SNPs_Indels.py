#!/usr/bin/env python3

# GWAS batch scripts for the 96-strain MAF15 SNP and InDel runs, all conditions.
# Writes, for each condition, <condition fastlmm dir>/codes/batch_scripts_with_double_id_snps_indels/
# (run_task_<VARIANT>.bash, submit_gwas_<VARIANT>.sh, submit_all_variants.sh, README.txt).
# Submission: through common_analysis_for_all/gwas_queue (queue.tsv lists every
# condition x variant; at most 4 variant arrays and 100 array tasks each at a time).
#
# SNPs   -> results/results_SNPs_96_MAF15_from_969_strains_with_double_id_final_filtered_for_ncount_less_than_5
#   full969Matrix...MAF5.NoSubTelo20kb.plink, 96 strains kept first (95 present,
#   SACE_YCY is not in the 969 file), then MAF >= 15% in those strains
#   (plink --keep-allele-order at both steps).
# Indels -> results/results_Indels_96_MAF15_with_double_id_final_filtered_for_ncount_less_than_5
#   indels.96samples (the InDel file used in Fig. 3, 96 strains, MAF >= 15%).
#   Not run for YPD.1 (its InDel results with this file already exist).
#
# --double-id IS used for both: these .fam files store strain names whole in both
# columns (e.g. SACE_GAL SACE_GAL). Without it runFaSTLMM splits phenotype names at
# '_' and all SACE_ strains are silently dropped (tested on YPD.1 YAL001C: 75 vs 82
# strains matched).
import math
import os

def create_batch_scripts_for_all_variants(variant_types, batch_dir, script_output_dir):
    """Create GWAS batch processing scripts for the given variant types"""

    # Slurm time limit per variant type = estimated batch time (genes per batch
    # x sec_per_gene x safety_factor, rounded up to whole hours) + buffer_hours.
    # On timeout, a task is resubmitted with extra_hours more, up to max_attempts.
    safety_factor = 2
    buffer_hours = 1
    min_hours = 4       # never request less than this
    extra_hours = 4
    max_attempts = 4
    # at most this many array tasks of one variant type run at the same time
    max_concurrent_tasks = 100

    
    # Create output directory
    os.makedirs(script_output_dir, exist_ok=True)
    
    
    print("CREATING GWAS BATCH SCRIPTS FOR ALL VARIANT TYPES")
    
    
    # Count batch files
    batch_files = [f for f in os.listdir(batch_dir) 
                   if f.startswith("phenotype_batch_") and f.endswith(".txt")
                   and "summary" not in f]
    num_batches = len(batch_files)
    genes_per_batch = max(
        sum(1 for line in open(os.path.join(batch_dir, f)) if line.strip())
        for f in batch_files
    ) if batch_files else 0

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

        est_hours = genes_per_batch * variant_config['sec_per_gene'] * safety_factor / 3600
        time_hours = max(math.ceil(est_hours) + buffer_hours, min_hours)
        print(f"Estimated time per batch: {est_hours:.2f} h "
              f"({genes_per_batch} genes x {variant_config['sec_per_gene']} s x {safety_factor}) "
              f"-> --time={time_hours}:00:00\n")

        # Create run-task.bash script
        run_task_script = os.path.join(script_output_dir, f"run_task_{variant_type}.bash")
        
        run_task_content = f"""#!/bin/bash
export PST_NUM_THREADS=1
# one math-library thread per task (the job has --cpus-per-task=1)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# Set variables
GENOTYPE="{genotype_path}"
OUTDIR="{output_dir}"
NPERM=100
THREADS=1

# Batch file directory
BATCH_DIR="{batch_dir}"
BATCH_FILE="${{BATCH_DIR}}/phenotype_batch_$(printf '%02d' $SLURM_ARRAY_TASK_ID).txt"

# Create output directory
mkdir -p $OUTDIR

# Run FaST-LMM on the node's local disk ($TMPDIR) instead of /projects: copy the
# genotype files there, write each gene's output/work files there, and copy only
# the finished gene folder to $OUTDIR. Slurm deletes $TMPDIR when the job ends.
LOCAL_DIR="${{TMPDIR:-/tmp}}/gwas_{variant_type}_${{SLURM_JOB_ID}}"
LOCAL_OUT="$LOCAL_DIR/out"
mkdir -p "$LOCAL_DIR/geno" "$LOCAL_OUT"
for ext in bed bim fam; do
    cp "$GENOTYPE.$ext" "$LOCAL_DIR/geno/"
done
LOCAL_GENOTYPE="$LOCAL_DIR/geno/$(basename "$GENOTYPE")"
echo "Local working directory: $LOCAL_DIR"

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
echo "Expected duration: ~$((BATCH_SIZE * {variant_config['sec_per_gene']} / 60)) minutes ({variant_config['sec_per_gene']} s per gene)"

# Read each phenotype file from the batch and process it
COUNTER=1
while IFS= read -r PHENOTYPE; do
    # Remove any trailing whitespace/newlines
    PHENOTYPE=$(echo "$PHENOTYPE" | tr -d '\\r\\n')
    
    # Skip empty lines
    if [ -z "$PHENOTYPE" ]; then
        continue
    fi
    
    GENE_NAME=$(basename "$PHENOTYPE" .norm.txt)
    
    echo
    echo "--- Processing file $COUNTER/$BATCH_SIZE: $GENE_NAME ---"
    echo "Phenotype file: $PHENOTYPE"
    
    # Check if file exists
    if [ ! -f "$PHENOTYPE" ]; then
        echo "ERROR: Phenotype file $PHENOTYPE does not exist!"
        COUNTER=$((COUNTER + 1))
        continue
    fi
    
    # runFaSTLMM names the output folder after the phenotype column header
    PHEN_NAME=$(head -1 "$PHENOTYPE" | cut -f2)
    PHEN_DIR="$OUTDIR/$PHEN_NAME"

    # On a resubmitted (timed-out) job: skip genes that already finished,
    # and remove the half-finished folder of the gene that was interrupted
    if [ -d "$PHEN_DIR" ]; then
        if grep -q "=== Done: $PHEN_NAME ===" "$PHEN_DIR/$PHEN_NAME.log" 2>/dev/null; then
            echo "SKIPPED: $GENE_NAME already completed"
            COUNTER=$((COUNTER + 1))
            continue
        fi
        echo "Removing incomplete output from an earlier attempt: $PHEN_DIR"
        rm -rf "$PHEN_DIR"
    fi
    # leftover of a copy to $OUTDIR that was interrupted
    rm -rf "$OUTDIR/.incoming_$PHEN_NAME"
    rm -rf "$LOCAL_OUT/$PHEN_NAME"

    # Run GWAS for this phenotype
    echo "Running GWAS for $GENE_NAME..."
    START_TIME=$(date +%s)

    python /scratch/aa9679/victor-material/Caudal2024_repro/runFaSTLMM.v4.0.py \\
        -g $LOCAL_GENOTYPE \\
        -p $PHENOTYPE \\
        -o $LOCAL_OUT \\
        --nb-permutations $NPERM \\
        --double-id \\
        --uncompressed \\
        --threads $THREADS

    EXIT_CODE=$?
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))

    if [ $EXIT_CODE -eq 0 ]; then
        # copy under a temporary name, then rename: a gene folder only appears
        # in $OUTDIR once it is complete
        cp -r "$LOCAL_OUT/$PHEN_NAME" "$OUTDIR/.incoming_$PHEN_NAME" && \\
            mv "$OUTDIR/.incoming_$PHEN_NAME" "$PHEN_DIR" && \\
            rm -rf "$LOCAL_OUT/$PHEN_NAME"
        echo "SUCCESS: $GENE_NAME completed in ${{DURATION}}s"
    else
        echo "FAILED: $GENE_NAME failed with exit code $EXIT_CODE after ${{DURATION}}s"
        # keep the failed gene's log for debugging
        rm -rf "$OUTDIR/.failed_$PHEN_NAME"
        [ -d "$LOCAL_OUT/$PHEN_NAME" ] && cp -r "$LOCAL_OUT/$PHEN_NAME" "$OUTDIR/.failed_$PHEN_NAME" && rm -rf "$LOCAL_OUT/$PHEN_NAME"
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
#SBATCH --time={time_hours}:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --account=torch_pr_130_general
#SBATCH --mail-type=FAIL,END
#SBATCH --array=1-{num_batches}%{max_concurrent_tasks}
#SBATCH --signal=B:USR1@600

# Timeout handling: Slurm sends USR1 10 min before the time limit. The trap
# resubmits only this array task with {extra_hours} more hours; the new run
# skips genes that already finished (see run_task_{variant_type}.bash).
TIME_HOURS=${{TIME_HOURS:-{time_hours}}}
ATTEMPT=${{ATTEMPT:-1}}
MAX_ATTEMPTS={max_attempts}

resubmit_on_timeout() {{
    echo "=== Array task $SLURM_ARRAY_TASK_ID (job $SLURM_JOB_ID) is about to reach its ${{TIME_HOURS}}h limit (attempt $ATTEMPT) ==="
    if [ "$ATTEMPT" -lt "$MAX_ATTEMPTS" ]; then
        NEW_HOURS=$((TIME_HOURS + {extra_hours}))
        sbatch --array=$SLURM_ARRAY_TASK_ID \\
               --time=${{NEW_HOURS}}:00:00 \\
               --export=ALL,TIME_HOURS=$NEW_HOURS,ATTEMPT=$((ATTEMPT + 1)) \\
               {slurm_script}
        echo "Resubmitted array task $SLURM_ARRAY_TASK_ID with --time=${{NEW_HOURS}}:00:00 (attempt $((ATTEMPT + 1)))"
    else
        echo "Reached MAX_ATTEMPTS=$MAX_ATTEMPTS, not resubmitting array task $SLURM_ARRAY_TASK_ID"
    fi
    kill -TERM "$TASK_PID" 2>/dev/null
    exit 1
}}
trap resubmit_on_timeout USR1

echo "Attempt $ATTEMPT, time limit ${{TIME_HOURS}}h"

# Load required modules
module purge

# Run the task script in the background so the USR1 trap can fire while it runs
/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/singularity/run-fastlmm-0.6.12.bash bash {run_task_script} &
TASK_PID=$!
wait "$TASK_PID"
"""
        
        with open(slurm_script, 'w') as f:
            f.write(slurm_content)
        
        os.chmod(slurm_script, 0o755)
        print(f"Created: {os.path.basename(slurm_script)}")
        
        job_info.append({
            'variant_type': variant_type,
            'batches': num_batches,
            'submit_script': slurm_script,
            'time_hours': time_hours,
            'est_hours': est_hours
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
        readme_content += f"\n  Output:   {config['outdir']}"
        readme_content += f"\n  Batches:  {info['batches']}"
        readme_content += (f"\n  Time:     --time={info['time_hours']}:00:00 "
                           f"(estimated {info['est_hours']:.2f} h per batch = "
                           f"{genes_per_batch} genes x {config['sec_per_gene']} s x {safety_factor})")

    submit_lines = "\n".join(f"   sbatch submit_gwas_{info['variant_type']}.sh"
                             for info in sorted(job_info, key=lambda x: x['variant_type']))
    readme_content += f"""

DIRECTORY STRUCTURE:
Each variant type has:
  - submit_gwas_[VARIANT].sh (SLURM submission script)
  - run_task_[VARIANT].bash (processing script)

USAGE:

1. Submit ALL variant types:
   bash submit_all_variants.sh

2. Submit individual variant type:
{submit_lines}

3. Monitor jobs:
   squeue -u $USER

4. Check results (one folder per gene inside each Output folder listed above):
   ls <Output folder>/<gene>/

Slurm logs (gwas_[VARIANT]_<jobid>_<task>.o/.e) are written to the directory
you submit from.

PARAMETERS:
- Permutations: 100
- Threads: 1 (--cpus-per-task=1)
- Memory: 10GB per job
- Time limit: per variant type, see above (minimum {min_hours} h)
- Array jobs: {num_batches} batches per variant type, up to {genes_per_batch} genes per batch,
  at most {max_concurrent_tasks} running at once (--array=1-{num_batches}%{max_concurrent_tasks})
- FaST-LMM runs on the node's local disk ($TMPDIR); only finished gene folders
  are copied to the Output folder (via a temporary .incoming_<gene> name)
- --double-id: used (FID = IID = full strain name, e.g. SACE_GAL)

TIMEOUT HANDLING:
- 10 min before the time limit Slurm sends USR1; the job resubmits only that
  array task with {extra_hours} more hours, up to {max_attempts} attempts.
- A rerun skips genes whose log contains "=== Done: <gene> ===", deletes the
  incomplete folder of the interrupted gene and reruns it.
- Resubmitting a whole array to the same Output folder also skips finished genes.
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

ROOT = "/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all"
BATCH_SCRIPTS_NAME = "batch_scripts_with_double_id_snps_indels"

SNP_GENOTYPE = '/projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/pan_transcriptome_paper/SNPs/SNPs_96_strain_filtering_no_change_name/SNPs_96_strains_MAF15_filtering/SNPs.96strains.maf15.plink'
INDEL_GENOTYPE = '/projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/vcf_files_sent_by_Anne/SNPs_indels/modified_snp_indel_plink_for_96strains/Indels_only_96/indels.96samples'

# variant type -> genotype, results folder name, sec_per_gene (wall time of one
# gene, 100 permutations, 1 thread; InDels scaled from 14 s for 1,033 variants)
RUNS = {
    'SNPs': {
        'sec_per_gene': 118,
        'genotype': SNP_GENOTYPE,
        'results': 'results_SNPs_96_MAF15_from_969_strains_with_double_id_final_filtered_for_ncount_less_than_5',
    },
    'Indels': {
        'sec_per_gene': 40,
        'genotype': INDEL_GENOTYPE,
        'results': 'results_Indels_96_MAF15_with_double_id_final_filtered_for_ncount_less_than_5',
    },
}

# condition -> (its fastlmm folder relative to ROOT, variant types to run)
CONDITIONS = {
    'YPD.1':              ('YPD.1/fastlmm',                               ['SNPs']),
    'YPD.2':              ('YPD.2/fastlmm',                               ['SNPs', 'Indels']),
    'Flu.1':              ('Flu.1/fastlmm',                               ['SNPs', 'Indels']),
    'Flu.2':              ('Flu.2/fastlmm',                               ['SNPs', 'Indels']),
    'merged_fluconazole': ('merged_fluconazole/new_fastlmm_march_26_2026', ['SNPs', 'Indels']),
    'sc2glu.2':           ('sc2glu.2/fastlmm',                            ['SNPs', 'Indels']),
    'scgal.1':            ('scgal.1/fastlmm',                             ['SNPs', 'Indels']),
    'merged_YPD':         ('merged_YPD/fastlmm',                          ['SNPs', 'Indels']),
}

if __name__ == "__main__":
    for condition, (fastlmm_dir, variants) in CONDITIONS.items():
        base = os.path.join(ROOT, fastlmm_dir)
        print(f"\n##### {condition}: {', '.join(variants)} #####")
        variant_types = {
            vt: {
                'sec_per_gene': RUNS[vt]['sec_per_gene'],
                'genotype': RUNS[vt]['genotype'],
                'outdir': os.path.join(base, "results", RUNS[vt]['results']),
            }
            for vt in variants
        }
        create_batch_scripts_for_all_variants(
            variant_types,
            batch_dir=os.path.join(base, "data", "pheno_batches"),
            script_output_dir=os.path.join(base, "codes", BATCH_SCRIPTS_NAME),
        )
