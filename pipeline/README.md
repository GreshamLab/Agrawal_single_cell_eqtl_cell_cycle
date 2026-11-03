# Upstream pipeline (before the figure code)

These scripts turn raw 10x reads into the inputs that the figure code reads: Cell Ranger counts,
demultiplexed strain identities (demuxlet), pseudobulk expression and FaST-LMM eQTL results.
They were run on the NYU HPC (Slurm) from one folder per dataset
(`common_analysis_for_all/<dataset>/`). The same script was used for every dataset with only the
dataset name and paths changed, so one copy is kept here; the table below lists what changes.

Datasets: YPD.1, YPD.2, Flu.1 (fluconazole 16 ug/ml), Flu.2 (8 ug/ml), sc2glu.2 (synthetic
complete with glucose), scgal.1 (synthetic complete with galactose).

## Order

| Step | Code | Output |
|---|---|---|
| 1. Align and count | `01_alignment/align.sh` (Cell Ranger, custom reference) | `alignment/<run>/run_count/outs/` |
| 2. Split BAM by cell barcode | `02_demultiplex/standard/split_cell_barcodes.py`, `split_bam.py`, `split_bam_main.sh` | one BAM per barcode batch |
| 3. Pileup | `02_demultiplex/standard/parallel_pileup.py`, `pileup_run.sh` (popscle dsc-pileup) | `pileup/` |
| 4. Demuxlet | `02_demultiplex/standard/demuxlet.py`, `run_demuxlet.sh` (popscle demuxlet); the merged table is written by the `find ... > merged_demuxlet_output_*.tsv` line in the launcher | `demuxlet/merged_demuxlet_output_*.tsv` |
| 5. Genotypes | `Figure.3/Fig.3_code_to_run_pseudobulk_and_fastlmm/PartA_Prepare_vcf_file/` (gene-CNVs without chr17: `make_no_chr17_plink.sh`) | PLINK files, 96 strains |
| 6. Pseudobulk | `Figure.3/.../PartB.../Step1.Run_Pseudobulk.Rmd` (strains with > 5 cells); sc2glu.2: `05_pseudobulk_sc2glu/0_Pseudobulk_fastlmm_sc2glu.2_only.Rmd` | `fastlmm/data/pheno_data/` |
| 7. FaST-LMM | `Figure.3/.../PartB.../Step2` (phenotype batches), `Step3` (SV / gene-CNV batch scripts), `Step3b` (SNP / InDel batch scripts), `runFaSTLMM.v4.0.py`; submitted with `03_fastlmm_launch/` | `fastlmm/results/` |
| 8. Gene-CNVs without chr17 | `04_cnv_no_chr17/` (batch scripts = gene-CNV scripts of step 7 with the no_chr17 genotype) | `results_CNVs_..._no_chr17` |
| 9. Significant eQTLs | `Figure.3/.../PartB.../Step5.Combine_signif_snps_output_thr_perm.py` | `all_signif_*` tables read by Figures 3-4 |

Other inputs made here:

- `05_pseudobulk_sc2glu/merged_outputs_from_sc2glu.1_and_sc2glu.2.Rmd`: pseudobulk of sc2glu.1 and
  sc2glu.2 together. Used by Figure 1C (synthetic complete panels) and the DESeq2 supplementary
  tables; run with `Figure.1/.../codes/0_pseudobulk_no_cell_filter.sbatch`.
- `06_cell_cycle_YPD.2/Fig5.cell-cycle-harmonic.Rmd`: harmonic fits for YPD.2, read by
  `Figure.5/.../Fig5.-part-II-compare-with-ypd.2.Rmd`. Cell-cycle times for YPD.2 come from the
  scripts in `Figure.5/.../part1-bash-code-to-assign-cell-cycle-time/` with the YPD.2 paths.

## What changes per dataset

| Dataset | Cell Ranger output | Demultiplexing | Pseudobulk |
|---|---|---|---|
| YPD.1 | `output_of_alignment_YPD1_with_clean_gtf` | standard | Step1 |
| YPD.2 | `output_of_alignment_YPD2_with_clean_gtf` | `02_demultiplex/YPD.2_parallel/` (subset-bam, parallel pileup) | Step1 |
| Flu.1 | `output_of_alignment_flu1_cust_ref` | standard | Step1 |
| Flu.2 | `output_of_alignment_flu2_cust_ref` | standard | Step1 |
| sc2glu.2 | `sc.2_with_cus_ref` | standard | `0_Pseudobulk_fastlmm_sc2glu.2_only.Rmd` |
| scgal.1 | `scgal_with_clean_gtf` | standard | Step1 |

FaST-LMM runs used in the paper (results folder names in `Step5`): SNPs 96_MAF15, InDels 96_MAF15
(YPD.1: the earlier `results_Indels_with_double_id_final_filtered_for_ncount_less_than_5` run on the
same genotype file), SVs from the new PLINK files, gene-CNVs without chr17.

## Notes

- Paths are absolute paths on the NYU HPC; change them to rerun elsewhere.
- Some scripts were edited after they ran, and only output folder names differ from the run:
  `02_demultiplex/standard/demuxlet.py` (YPD.1 copy) writes `demuxlet_YPD1_out_1SD_default` (the run used `noSD`);
  the YPD.2 `split_bam_new_code.py` writes `bams_Flu.1_new` (the run wrote `bams_YPD.2_new`).
- `runFaSTLMM.v4.0.py` runs inside the FaST-LMM 0.6.12 Singularity image
  (`03_fastlmm_launch/run-fastlmm-0.6.12.bash`); the image is not included.
