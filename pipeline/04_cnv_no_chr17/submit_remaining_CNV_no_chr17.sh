#!/bin/bash
#SBATCH --job-name=submit_CNV_noChr17
#SBATCH --output=/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/submit_remaining_CNV_no_chr17_%j.log
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --account=torch_pr_130_general
# Submits the CNV FaST-LMM re-runs without chr17 for the remaining conditions, one array at a
# time, once the user's queue (limit 3000 array tasks) has room for the whole array + 50 spare
# slots for timeout resubmits.
cd /projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all
B=batch_scripts_with_double_id_is_not_needed_final_filtered_for_ncount_less_than_5_run_from_new_plink_files_no_chr17
for R in merged_fluconazole/new_fastlmm_march_26_2026 sc2glu.1/fastlmm_sc2.glu.1_only sc2glu.2/fastlmm sc2glu.2/fastlmm_sc2.glu.2_only scgal.1/fastlmm delta_expression_eQTL/Flu16_minus_YPDint/fastlmm delta_expression_eQTL/Flu8_minus_YPDint/fastlmm delta_expression_eQTL/YPD1_minus_YPD2/fastlmm; do
  D=$R/codes/$B
  size=$(grep -o 'array=1-[0-9]*' $D/submit_gwas_CNVs.sh | cut -d- -f2)
  while true; do
    n=$(squeue -u $USER -h -r 2>/dev/null | wc -l)
    if [ $((n + size + 50)) -le 3000 ]; then
      out=$(cd $D && sbatch submit_gwas_CNVs.sh 2>&1)
      case "$out" in *Submitted*) echo "$(date '+%F %H:%M') $R: $out"; break;; esac
      echo "$(date '+%F %H:%M') $R: retry ($out)"
    fi
    sleep 120
  done
done
echo "ALL SUBMITTED $(date '+%F %H:%M')"
