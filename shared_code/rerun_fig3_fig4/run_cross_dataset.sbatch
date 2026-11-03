#!/bin/bash
#SBATCH --account=torch_pr_130_general
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --output=logs/%x_%j.out

# Analyses that combine datasets; run after every per-dataset job has finished.
P=/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle
R=$P/shared_code/rerun_fig3_fig4/render_rmd.sh
F3S=$P/Figure.3/Fig.3-supplementary/reproducibility_of_eQTL_mapping_between_conditions/codes
F4=$P/Figure.4/Fig.4-PartII-compare-conditions/codes
module purge
module load r/4.6.1
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

# Figure 3 supplementary: YPD.1 vs YPD.2
bash $R $F3S/Fig3-supp-YPD1-YPD2-reproducibility.Rmd
bash $R $F3S/Fig.3-supp-YPD.1-YPD.2-expression-variation.Rmd

# Figure 4 Part II: gene-level min P (GO universe), then the comparisons
(cd $F4 && echo "$(date '+%F %T')  START  extract_gene_level_minP.py" &&
 python3 extract_gene_level_minP.py > extract_gene_level_minP.render_2026-11-01_arial.log 2>&1 &&
 echo "$(date '+%F %T')  DONE  extract_gene_level_minP.py") || echo "FAILED extract_gene_level_minP.py"
for rmd in Fig4.PartII-compare-conditions-replicates.Rmd \
           Fig4.effect-size-comparison-replicates.Rmd \
           Fig4.effect-size-comparison-flu-doses-vs-baselines.Rmd \
           Fig4.flu-eQTLs-in-baselines-concordance-and-GO.Rmd \
           Fig4.mashr-condition-specific-eQTLs.Rmd \
           Fig4.mashr-per-variant-type.Rmd \
           Fig4.condition-induced-genes-eQTLs.Rmd; do
  bash $R $F4/$rmd
done

# Figure 4 supplementary: single-cell vs bulk eQTL effect sizes (sc2glu.2)
bash $R "$P/Figure.4/Fig.4-supplementary/Agreement_of_sc2%glucose_10x_vs_bulk/Step3.effect-size-correlation-bulk-vs-sc2glu.2-only.Rmd"
echo "cross-dataset job finished $(date)"
