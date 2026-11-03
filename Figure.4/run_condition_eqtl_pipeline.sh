#!/bin/bash
# Per-condition eQTL pipeline (same design as Figure.3 / YPD.1):
#   Step1 classification Rmd: per-gene LD clumping in R (r2 > 0.6 within 50 kb,
#         shared_code/ld_clump_and_local.R), annotation, cis/trans (2 kb) and local/distant
#         (25 kb; CNVs same ORF), panels B-D
#   Step2 (panel A: cis/trans map and local/distant map), Step3 (single-cell verification of local eQTLs)
# 2026-10-16: PLINK clumping removed. It clumped all genes together with a 2-kb window; the
#   per-gene clumping now runs inside Step1 (cached in results/<condition>/clumped/).
# Usage: bash run_condition_eqtl_pipeline.sh <condition>   e.g. Flu.1, Flu_integrated
#        bash run_condition_eqtl_pipeline.sh all           (all conditions in parallel)
set -u
F4=/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle/Figure.4
V=/projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs
CONDITIONS="Flu.1 Flu.2 Flu_integrated sc2glu_integrated sc2glu.1_only sc2glu.2_only scgal.1 YPD.2 YPD_integrated"

module purge
module load r/4.6.1

render() { Rscript -e "rmarkdown::render('$1')"; }

run_one() {
  local d=$1 root=$F4/Fig.4-$1
  local codes=$root/codes logs=$root/results/logs
  mkdir -p "$logs"
  cd "$codes" || return 1

  # knit the three Rmds (Step1 does the clumping)
  render Step1_*.Rmd > "$logs/step1.log" 2>&1 || { echo "$d: Step1 FAILED"; return 1; }
  render Step2_*.Rmd > "$logs/step2.log" 2>&1 || { echo "$d: Step2 FAILED"; return 1; }
  render Step3_*.Rmd > "$logs/step3.log" 2>&1 || { echo "$d: Step3 FAILED"; return 1; }
  echo "$d: DONE"
}

if [ "${1:-}" = "all" ]; then
  for d in $CONDITIONS; do run_one "$d" & done; wait
else
  run_one "$1"
fi
