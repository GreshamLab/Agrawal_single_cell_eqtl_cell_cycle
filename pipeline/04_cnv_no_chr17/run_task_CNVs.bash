#!/bin/bash
export PST_NUM_THREADS=1

# Set variables
GENOTYPE="/projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/structural_variant_vcf/SVcalls2/SVcalls/CNVs_96_strains/CNV_96strains_maf15/no_chr17/CNVs.96strains.maf15.no_chr17.plink"
OUTDIR="/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/YPD.1/fastlmm/results/results_CNVs_with_double_id_is_not_needed_final_filtered_for_ncount_less_than_5_run_from_new_plink_files_no_chr17"
NPERM=100
THREADS=1

# Batch file directory
BATCH_DIR="/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/YPD.1/fastlmm/data/pheno_batches"
BATCH_FILE="${BATCH_DIR}/phenotype_batch_$(printf '%02d' $SLURM_ARRAY_TASK_ID).txt"

# Create output directory
mkdir -p $OUTDIR

echo "=== GWAS BATCH JOB $SLURM_ARRAY_TASK_ID FOR CNVs ===";
echo "Variant type: CNVs"
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
echo "Expected duration: ~$((BATCH_SIZE * 16 / 60)) minutes (16 s per gene)"

# Read each phenotype file from the batch and process it
COUNTER=1
while IFS= read -r PHENOTYPE; do
    # Remove any trailing whitespace/newlines
    PHENOTYPE=$(echo "$PHENOTYPE" | tr -d '\r\n')
    
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

    # Run GWAS for this phenotype
    echo "Running GWAS for $GENE_NAME..."
    START_TIME=$(date +%s)
    
    python /scratch/aa9679/victor-material/Caudal2024_repro/runFaSTLMM.v4.0.py \
        -g $GENOTYPE \
        -p $PHENOTYPE \
        -o $OUTDIR \
        --nb-permutations $NPERM \
        --uncompressed \
        --threads $THREADS
    
    EXIT_CODE=$?
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    
    if [ $EXIT_CODE -eq 0 ]; then
        echo "SUCCESS: $GENE_NAME completed in ${DURATION}s"
    else
        echo "FAILED: $GENE_NAME failed with exit code $EXIT_CODE after ${DURATION}s"
    fi
    
    COUNTER=$((COUNTER + 1))
    
done < "$BATCH_FILE"

echo
echo "=== BATCH $SLURM_ARRAY_TASK_ID COMPLETED FOR CNVs ==="
echo "Finished at: $(date)"
echo "Processed $((COUNTER - 1)) phenotype files"

# Create a completion marker
touch "$OUTDIR/batch_${SLURM_ARRAY_TASK_ID}_completed.marker"
