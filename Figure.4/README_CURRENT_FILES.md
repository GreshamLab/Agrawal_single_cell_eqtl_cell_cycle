# Figure 4: which files are current (updated 2026-10-09)

Files named `OLD_*`, and anything in an `OLD_2026-09-30_run/` folder, are outdated. They use the
old integrated datasets or old variant files (full969 SNPs, CNV_969strains, SV.1087), or are
outputs of the 30 Sept run. Do not use them for the figure.

All datasets below use the YPD.1 variant files (SNPs.96strains.maf15, indels.96samples,
CNVs.96strains.maf15, SV.96strains.maf15) and the YPD.1 pseudobulk steps.

## Update 2026-10-16: per-gene LD clumping and local/distant eQTLs

- Clumping now runs inside Step1 in R (`../shared_code/ld_clump_and_local.R`), **separately for
  every gene**: significant variants with r2 > 0.6 within 50 kb are one signal. The earlier PLINK
  step (`--clump-kb 2`) clumped all genes together (its input had no gene column) and is removed
  from `run_condition_eqtl_pipeline.sh`. Window chosen from the LD decay of the 96 strains
  (`Figure.3/.../codes/LD_decay_96_strains.Rmd`).
- eQTLs are classified as **local** (any variant of the clump within 25 kb of the gene; SV span;
  translocations either breakpoint; gene-based CNVs only the same ORF) or **distant**. Figures
  in Step1 use local/distant; Step2 draws two panel A maps (`*_panelA_cis_trans_2x2` and
  `*_panelA_local_distant_2x2`); Step3 verifies local eQTLs (`PanelB_verify_<type>_local.csv`).
- The annotated tables keep their names (`*_annotated_cis_trans.txt`), keep the old `cis_trans`
  column (lead variant within 2 kb) and add `local_distant`, `local_distance_bp`,
  `clump_n_members`, `clump_members`.
- Older outputs (PLINK `clumped_results_*.clumped`, `Fig1A_effect_size_cis_trans.tiff`,
  `PanelB_verify_*_cis.csv`) come from the previous clumping and are outdated.
- The Part II comparison Rmds read these annotation tables, so they change when re-knitted.

## 1. Individual figures per condition

| Code | Role |
|---|---|
| `conditions.tsv` | condition -> FaST-LMM results, gene_location, Seurat object |
| `generate_condition_rmds.py` | writes the three Rmds per condition from the YPD.1 Figure 3 Rmds |
| `run_condition_eqtl_pipeline.sh` | Step1 -> PLINK clumping -> Step1/2/3 knits for one condition |
| `run_condition_eqtl_pipeline.sbatch` | Slurm wrapper: `sbatch --job-name=<c> run_condition_eqtl_pipeline.sbatch <c>` |
| `Fig.4-<condition>/codes/Step1_*.Rmd` | eQTL classification, panels B-D |
| `Fig.4-<condition>/codes/Step2_*.Rmd` | panel A, cis/trans eQTL map |
| `Fig.4-<condition>/codes/Step3_*.Rmd` | raw single-cell verification of cis-eQTLs |

Conditions: YPD.2, YPD_integrated, Flu.1, Flu.2, Flu_integrated, sc2glu.1_only,
sc2glu.2_only, sc2glu_integrated, scgal.1 (YPD.1 is Figure 3).

Outputs per condition:
- eQTL map: `Fig.4-<c>/results/<c>/eqtl_plot/figures/Fig4_<c>_panelA_cis_trans_2x2.{tiff,png}`
- panels B-D: `Fig.4-<c>/results/<c>/eqtl_plot/figures/`
- raw single-cell plots: `Fig.4-<c>/results/raw_single_cell_plots/`
- annotated eQTL tables: `Fig.4-<c>/results/<c>/annotations/`

## 2. Comparing conditions (`Fig.4-PartII-compare-conditions/codes/`)

| Code | Datasets | Output (`../results/`) |
|---|---|---|
| `Fig4.PartII-compare-conditions-best-replicates.Rmd` | YPD.1, sc2glu.2, Flu.1, Flu.2 | `best_replicates/` |
| `Fig4.effect-size-comparison-best-replicates.Rmd` | YPD.1, sc2glu.2, Flu.1, Flu.2 (all ordered pairs) | `best_replicates_concordance/` |
| `Fig4.effect-size-comparison-flu-doses-vs-baselines.Rmd` | Flu.1, Flu.2 vs YPD.1, YPD.2, sc2glu.2, per variant type | `effect_size_flu_doses_vs_baselines/` |
| `Fig4.effect-size-comparison-replicates.Rmd` | all 6 datasets (input for the Rmd above) | `replicates_concordance/` |
| `Fig4.PartII-compare-conditions-replicates.Rmd` | all 6 datasets, 5 dataset sets | `replicates/` |
| `render_best_replicates.sbatch` | knits any of these: `sbatch render_best_replicates.sbatch <Rmd>` | |

In `Fig.4-PartII-compare-conditions/results/`, all `OLD_*` folders come from the OLD
integrated-dataset Rmds (`OLD_integrated_top_level/` holds their loose TIFF/CSV/RDS files).
