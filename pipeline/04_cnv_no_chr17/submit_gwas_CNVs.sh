#!/bin/bash
#
#SBATCH --verbose
#SBATCH --job-name=GWAS_CNVs_noChr17
#SBATCH --output=gwas_CNVs_%A_%a.o
#SBATCH --error=gwas_CNVs_%A_%a.e
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --account=torch_pr_130_general
#SBATCH --array=1-531
#SBATCH --signal=B:USR1@600

# Timeout handling: Slurm sends USR1 10 min before the time limit. The trap
# resubmits only this array task with 4 more hours; the new run
# skips genes that already finished (see run_task_CNVs.bash).
TIME_HOURS=${TIME_HOURS:-4}
ATTEMPT=${ATTEMPT:-1}
MAX_ATTEMPTS=4

resubmit_on_timeout() {
    echo "=== Array task $SLURM_ARRAY_TASK_ID (job $SLURM_JOB_ID) is about to reach its ${TIME_HOURS}h limit (attempt $ATTEMPT) ==="
    if [ "$ATTEMPT" -lt "$MAX_ATTEMPTS" ]; then
        NEW_HOURS=$((TIME_HOURS + 4))
        sbatch --array=$SLURM_ARRAY_TASK_ID \
               --time=${NEW_HOURS}:00:00 \
               --export=ALL,TIME_HOURS=$NEW_HOURS,ATTEMPT=$((ATTEMPT + 1)) \
               /projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/YPD.1/fastlmm/codes/batch_scripts_with_double_id_is_not_needed_final_filtered_for_ncount_less_than_5_run_from_new_plink_files_no_chr17/submit_gwas_CNVs.sh
        echo "Resubmitted array task $SLURM_ARRAY_TASK_ID with --time=${NEW_HOURS}:00:00 (attempt $((ATTEMPT + 1)))"
    else
        echo "Reached MAX_ATTEMPTS=$MAX_ATTEMPTS, not resubmitting array task $SLURM_ARRAY_TASK_ID"
    fi
    kill -TERM "$TASK_PID" 2>/dev/null
    exit 1
}
trap resubmit_on_timeout USR1

echo "Attempt $ATTEMPT, time limit ${TIME_HOURS}h"

# Load required modules
module purge

# Run the task script in the background so the USR1 trap can fire while it runs
/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/singularity/run-fastlmm-0.6.12.bash bash /projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/YPD.1/fastlmm/codes/batch_scripts_with_double_id_is_not_needed_final_filtered_for_ncount_less_than_5_run_from_new_plink_files_no_chr17/run_task_CNVs.bash &
TASK_PID=$!
wait "$TASK_PID"
