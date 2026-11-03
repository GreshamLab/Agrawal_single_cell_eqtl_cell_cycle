#!/bin/bash
#SBATCH --job-name=render_Fig4_noChr17
#SBATCH --output=/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle/Figure.4/render_all_Fig4_no_chr17_%j.log
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --account=torch_pr_130_general
#SBATCH --mail-type=FAIL,END
# Re-renders Figure 4 after the CNV FaST-LMM re-run without chr17 (2026-10-23):
# every condition Step1 -> Step2 -> Step3, then the gene-level minP table, then PartII.
# A failed Rmd is reported and the rest continue; Step2/Step3 of a condition are skipped
# when its Step1 fails (they read Step1's annotation tables).

F4=/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle/Figure.4
R="singularity exec --overlay /share/apps/images/r-4.6.1-ubuntu-26.04.sqf:ro /share/apps/images/ubuntu-26.04.sif /ext3/r/4.6.1/bin/Rscript"
FAILED=()

render() {   # $1 = directory, $2 = Rmd
    local log="${2%.Rmd}.render_no_chr17.log"
    echo "$(date '+%F %H:%M') START $1/$2"
    (cd "$1" && $R -e "rmarkdown::render('$2')" > "$log" 2>&1)
    local st=$?
    if [ $st -eq 0 ]; then echo "$(date '+%F %H:%M') OK    $2"; else echo "$(date '+%F %H:%M') FAIL  $2 (see $1/$log)"; FAILED+=("$2"); fi
    return $st
}

for C in YPD.2 YPD_integrated Flu.1 Flu.2 Flu_integrated sc2glu.1_only sc2glu.2_only sc2glu_integrated scgal.1; do
    D=$F4/Fig.4-$C/codes
    if render $D Step1_Fig4_panel_B_C_D_variant_eqtl_classification_$C.Rmd; then
        render $D Step2_Fig4-panel-A-cis_trans_eqtl_mapping_$C.Rmd
        render $D Step3_Fig4-Panel-B-${C}_raw_eqtl_verification.Rmd
    else
        echo "  skipping Step2/Step3 of $C"
    fi
done

P=$F4/Fig.4-PartII-compare-conditions/codes
echo "$(date '+%F %H:%M') START extract_gene_level_minP.py"
(cd $P && python3 extract_gene_level_minP.py > extract_gene_level_minP.no_chr17.log 2>&1) \
    && echo "$(date '+%F %H:%M') OK    extract_gene_level_minP.py" \
    || { echo "$(date '+%F %H:%M') FAIL  extract_gene_level_minP.py"; FAILED+=("extract_gene_level_minP.py"); }

for f in Fig4.PartII-compare-conditions-replicates.Rmd Fig4.PartII-compare-conditions-best-replicates.Rmd \
         Fig4.effect-size-comparison-replicates.Rmd Fig4.effect-size-comparison-flu-doses-vs-baselines.Rmd \
         Fig4.flu-eQTLs-in-baselines-concordance-and-GO.Rmd Fig4.mashr-condition-specific-eQTLs.Rmd; do
    render $P $f
done

echo
echo "$(date '+%F %H:%M') DONE. Failed: ${#FAILED[@]} ${FAILED[*]}"
