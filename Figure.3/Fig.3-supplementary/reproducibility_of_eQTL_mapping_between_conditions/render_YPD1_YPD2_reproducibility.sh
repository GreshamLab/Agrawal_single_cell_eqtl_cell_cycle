#!/bin/bash
#SBATCH --job-name=render_repro_YPD1_YPD2
#SBATCH --output=/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle/Figure.3/Fig.3-supplementary/reproducibility_of_eQTL_mapping_between_conditions/render_YPD1_YPD2_reproducibility_%j.log
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=96gb
#SBATCH --account=torch_pr_130_general
#SBATCH --mail-type=FAIL,END
# eQTL reproducibility first: it writes shared_eQTLs_YPD1_vs_YPD2_ALL_variants.csv, which the
# expression-variation Rmd reads.
set -e
cd /projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle/Figure.3/Fig.3-supplementary/reproducibility_of_eQTL_mapping_between_conditions/codes
R="singularity exec --overlay /share/apps/images/r-4.6.1-ubuntu-26.04.sqf:ro /share/apps/images/ubuntu-26.04.sif /ext3/r/4.6.1/bin/Rscript"
for f in Fig3-supp-YPD1-YPD2-reproducibility.Rmd Fig.3-supp-YPD.1-YPD.2-expression-variation.Rmd; do
  echo "$(date '+%F %H:%M') START $f"
  $R -e "rmarkdown::render('$f')" > "${f%.Rmd}.render.log" 2>&1
  echo "$(date '+%F %H:%M') OK    $f"
done
echo "$(date '+%F %H:%M') DONE"
