# ============================================================
# Collect every paper-relevant result table (Figures 1-5 + supplements)
# into ONE Excel workbook, one sheet per table, with a README index sheet.
#
# Run AFTER 01_export_DESeq2_tables.R (needs the DESeq2 tables it writes).
# Output: Supplementary_Tables/results/Supplementary_Tables_all_figures.xlsx
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
})

PROJ <- "/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle"
OUT  <- file.path(PROJ, "Supplementary_Tables/results/Supplementary_Tables_all_figures.xlsx")

F1   <- "Figure.1"
F2M  <- "Figure.2/Fig.2-Panel-A-B-C-D-strain-clade-gene-expression-comparison/results"
F2S  <- "Figure.2/Fig.2_supplementary/results"
F3M  <- "Figure.3/Fig.3-Panel-A-C-D-cis-and-trans-eqtl-YPD-variant-type/results"
F3S  <- "Figure.3/Fig.3-supplementary"
F4M  <- "Figure.4/Fig.4-PartII-compare-conditions/results"
F4S  <- "Figure.4/Fig.4-supplementary/results"
F5   <- "Figure.5/2.cell_cycle_time_assignment_to_single_cells/results"
DE   <- "Supplementary_Tables/results/DESeq2_tables"
DE_CONDS <- c("YPD.1", "YPD.2", "Flu.1", "Flu.2", "sc2glu.2", "scgal.1")

rd <- function(path) {
  full <- file.path(PROJ, path)
  if (grepl("\\.csv$", path)) read_csv(full, show_col_types = FALSE, guess_max = 1e5)
  else read_tsv(full, show_col_types = FALSE, guess_max = 1e5)
}
# Read several files and stack them with an identifying column
stack <- function(paths, ids, id_col = "condition") {
  map2(paths, ids, ~ rd(.x) %>% mutate(!!id_col := .y, .before = 1)) %>% bind_rows()
}

# ── Sheet registry: name (<=31 chars), figure, description, data, sources ──
sheets <- list()
add <- function(name, figure, description, data, sources) {
  stopifnot(nchar(name) <= 31, !name %in% names(sheets))
  sheets[[name]] <<- list(figure = figure, description = description,
                          data = data, sources = paste(sources, collapse = "; "))
}

# ---------------- Figure 1 ----------------
add("F1_strains", "Fig. 1A / Table",
    "96 wild isolates used for multiplexed scRNA-seq: origin, ploidy, zygosity, clade (1011 Yeast Genomes).",
    rd(file.path(F1, "Fig.1-all-experiment-table/StrainsforSC-RNAseq.csv")),
    file.path(F1, "Fig.1-all-experiment-table/StrainsforSC-RNAseq.csv"))
f <- file.path(F1, "Fig.1-all-experiment-table/results/Table1_variant_type_subtype_summary.csv")
add("F1_variant_types", "Table 1", "Number of variants by type and subtype across the strain panel.", rd(f), f)
f <- file.path(F1, "Fig.1-all-experiment-table/results/Table2_pairwise_diversity_summary.csv")
add("F1_pairwise_diversity", "Table 2", "Pairwise genetic diversity between strains, per variant type.", rd(f), f)
f <- file.path(F1, "Fig.1-all-experiment-table/results/Supplementary_Table_QC_summary.csv")
add("F1_scRNAseq_QC", "Supp. Table", "scRNA-seq run QC per condition/replicate: cells, reads, genes, demuxlet singlet/doublet rates.", rd(f), f)
f <- file.path(F1, "Fig.1-Panel-B-demuxlet-accurately-detects-S288C/results/Fig1B_source_data.csv")
add("F1B_demuxlet_accuracy", "Fig. 1B", "Demuxlet recall/precision for S288C (DGY3) cells validated by KANMX/NATMX marker reads.", rd(f), f)
f <- file.path(F1, "Fig.1-Panel-C-scRNAseq-reproducibility_agreement/results/Fig1C_summary_stats.csv")
add("F1C_replicate_concordance", "Fig. 1C", "Spearman correlation of pseudobulk expression between biological replicates.", rd(f), f)
f <- file.path(F1, "Fig.1-Panel-C-scRNAseq-reproducibility_agreement/results/Fig1C_negative_control_summary_stats.csv")
add("F1C_negative_control", "Fig. 1C", "Negative-control comparisons (across conditions) for replicate concordance.", rd(f), f)
f <- file.path(F1, "Fig.1-Panel-C-scRNAseq-reproducibility_agreement/results/Fig1D_strain_concordance_summary_stats.csv")
add("F1D_strain_concordance", "Fig. 1D", "Same-strain vs different-strain expression correlation between replicates (Wilcoxon test).", rd(f), f)
f <- file.path(F1, "Fig.1_supplementary/results/FigS4-5_per_strain.csv")
add("F1S_per_strain_repro", "Supp. Fig. 4-5", "Per-strain replicate correlation with cell numbers and UMI depth.", rd(f), f)
f <- file.path(F1, "Fig.1_supplementary/results/FigS4-6_summary_stats.csv")
add("F1S_repro_summary", "Supp. Fig. 4-6", "Summary of per-strain and per-gene replicate reproducibility.", rd(f), f)
f <- file.path(F1, "Fig.1_supplementary/results/FigS6_per_gene_correlations.csv")
add("F1S_per_gene_corr", "Supp. Fig. 6", "Per-gene replicate correlation across strains.", rd(f), f)
f <- file.path(F1, "Fig.1_supplementary/results/SuppFig2_run_QC_summary.csv")
add("F1S_run_QC", "Supp. Fig. 2", "Per-run median reads, SNPs per cell and total cells.", rd(f), f)
f <- file.path(F1, "Fig.1_supplementary/results/SuppFig2_reproducibility_stats.csv")
add("F1S_genotype_abund_repro", "Supp. Fig. 2", "Replicate reproducibility of strain (genotype) abundance.", rd(f), f)
f <- file.path(F1, "Fig.1_supplementary/results/SuppFig3_within_clade_summary.csv")
add("F1S_ambiguous_in_clade", "Supp. Fig. 3", "Ambiguous demuxlet calls whose candidate strains fall within one clade.", rd(f), f)

# ---------------- Figure 2 ----------------
de_all <- map(DE_CONDS, ~ rd(file.path(DE, paste0("DESeq2_clade_vs_rest_all_genes_", .x, ".csv")))) %>% bind_rows()
add("F2_DESeq2_YPD1_all_genes", "Fig. 2A-C",
    "DESeq2 on pseudobulk counts, YPD.1: each clade (>=2 strains) vs all other strains, design ~group. All tested genes. DEG = padj<0.05 & |log2FC|>0.3.",
    de_all %>% filter(condition == "YPD.1"),
    file.path(DE, "DESeq2_clade_vs_rest_all_genes_YPD.1.csv"))
add("F2_DEGs_sig_all_conditions", "Fig. 2 / Supp. Fig. 2",
    "Significant DEGs only (padj<0.05 & |log2FC|>0.3) for every condition and clade. Full per-condition tables are in the DESeq2_tables folder.",
    de_all %>% filter(direction != "NS"),
    file.path(DE, paste0("DESeq2_clade_vs_rest_all_genes_", DE_CONDS, ".csv")))
add("F2_DEG_counts_per_clade", "Fig. 2 / Supp. Fig. 2", "Number of up/down-regulated genes per clade, per condition.",
    map(file.path(F2S, DE_CONDS, "ora_tables", paste0("clade_summary_", DE_CONDS, ".csv")), rd) %>%
      bind_rows(),   # files already carry a condition column
    file.path(F2S, "<condition>/ora_tables/clade_summary_<condition>.csv"))
f <- file.path(F2M, "Fig2_clade_GO_ORA_all.csv")
add("F2_GO_ORA_YPD1", "Fig. 2 GO dotplots", "GO over-representation (clusterProfiler) of up/down DEGs per clade, YPD.1.", rd(f), f)
f <- file.path(F2M, "Fig2_clade_GSEA_YPD1.csv")
add("F2_GSEA_YPD1", "Fig. 2", "GO GSEA on DESeq2 ranked genes per clade, YPD.1.", rd(f), f)
oc <- setdiff(DE_CONDS, "YPD.1")
add("F2S_GO_ORA_other_conds", "Supp. Fig. 2", "GO over-representation of DEGs per clade for YPD.2, Flu.1, Flu.2, sc2glu.2, scgal.1.",
    stack(file.path(F2S, oc, "ora_tables", paste0("GO_ORA_", oc, ".csv")), oc),
    file.path(F2S, "<condition>/ora_tables/GO_ORA_<condition>.csv"))
f <- file.path(F2S, "FigS_Fig2_strains_genes_conditions_summary.csv")
add("F2S_condition_summary", "Supp. Fig. 2", "Per condition: strains, genes detected, clades tested, total DEGs.", rd(f), f)
f <- file.path(F2S, "FigS_Fig2_regulon_selection_log.csv")
add("F2S_regulon_genes", "Fig. 2C / Supp.", "How GAL / NCR / stress regulon genes were selected for the focused panels.", rd(f), f)

# ---------------- Figure 3 (YPD.1 eQTL) + per-condition eQTL summaries ----------------
f <- file.path(F3M, "YPD.1/eqtl_plot/annotations/cis_trans_eQTL_mapping_all_variants_combined.txt")
add("F3_eQTL_YPD1_cis_trans", "Fig. 3A", "All clumped YPD.1 eQTLs (SNPs, InDels, CNVs, SVs) with target gene, causal gene, cis/trans call, effect size, P.", rd(f), f)

# Per-dataset eQTL outputs: YPD.1 is Figure 3, the other five paper datasets are
# Figure.4/Fig.4-<dataset>/ (sc2glu.2 = the sc2glu.2_only run; sc2glu.1 and the
# integrated datasets are not in the paper tables)
F4_DATASETS <- c(YPD.2 = "YPD.2", Flu.1 = "Flu.1", Flu.2 = "Flu.2",
                 sc2glu.2 = "sc2glu.2_only", scgal.1 = "scgal.1")
cond_sum <- c(YPD.1 = file.path(F3M, "YPD.1/eqtl_plot/data_summary/%s.txt"),
              setNames(file.path("Figure.4", paste0("Fig.4-", F4_DATASETS), "results", F4_DATASETS,
                                 "eqtl_plot/data_summary/%s.txt"), names(F4_DATASETS)))
sum_sheet <- function(stem) {
  p <- sprintf(cond_sum, stem); stack(p, names(cond_sum))
}
# Stack files whose column types differ (e.g. SNP vs SV tables): read as text, then re-guess types
stack_tc <- function(paths, ids, id_col = "condition") {
  map2(paths, ids, ~ read_tsv(file.path(PROJ, .x), col_types = cols(.default = "c")) %>%
         mutate(!!id_col := .y, .before = 1)) %>%
    bind_rows() %>% type_convert(guess_integer = FALSE)
}
# Part II tables cover sc2glu.1 / integrated datasets too; keep only the 6 paper datasets
PAPER_DATASETS <- c("YPD.1", "YPD.2", "Flu.1_16ug", "Flu.2_8ug", "sc2glu.2", "scgal.1")
paper_pairs <- function(d) d %>% filter(Condition1 %in% PAPER_DATASETS, Condition2 %in% PAPER_DATASETS)
add("F3_eQTL_variant_summary", "Fig. 3C / 4", "eQTL counts and percentages by variant type, per condition.",
    sum_sheet("variant_summary_stats"), sprintf(cond_sum, "variant_summary_stats"))
add("F3_eQTL_odds_ratio", "Fig. 3C / 4", "Fisher odds-ratio enrichment of each variant type among eQTLs, per condition.",
    sum_sheet("odds_ratio_enrichment"), sprintf(cond_sum, "odds_ratio_enrichment"))
add("F3_eQTL_pleiotropy", "Fig. 3D / 4", "Pleiotropy: number of genes affected per eQTL variant, by variant type, per condition.",
    sum_sheet("pleiotropy_summary"), sprintf(cond_sum, "pleiotropy_summary"))
add("F3_eQTL_pleiotropy_subtype", "Fig. 3D / 4", "Pleiotropy by variant subtype (e.g. DEL/INS, DUP), per condition.",
    sum_sheet("pleiotropy_by_subtype_summary"), sprintf(cond_sum, "pleiotropy_by_subtype_summary"))

f <- file.path(F3M, "raw_single_cell_plots/PanelB_validation_summary_all_types.csv")
add("F3B_local_validation_summary", "Fig. 3B", "Single-cell raw-data validation of local eQTLs: success rate and effect-direction agreement with FaST-LMM.", rd(f), f)
vt <- c("SNPs", "InDels", "CNVs", "SVs")
p <- file.path(F3M, "raw_single_cell_plots", paste0("PanelB_verify_", vt, "_local.csv"))
add("F3B_local_validation_detail", "Fig. 3B", "Per local eQTL single-cell verification (beta, P) vs FaST-LMM beta/P.",
    stack(p, vt, "variant_type"), p)

f <- file.path(F3S, "reproducibility_of_eQTL_mapping_between_conditions/results/Detection_Rate_Summary.csv")
add("F3S_YPD1v2_detection", "Supp. Fig. 3", "Fraction of YPD.1 eQTLs found in YPD.2 lookup, by variant type.", rd(f), f)
f <- file.path(F3S, "reproducibility_of_eQTL_mapping_between_conditions/results/Effect_Size_Concordance_Summary.csv")
add("F3S_YPD1v2_effect_conc", "Supp. Fig. 3", "Effect-size sign concordance and correlation, YPD.1 eQTLs looked up in YPD.2.", rd(f), f)
f <- file.path(F3S, "reproducibility_of_eQTL_mapping_between_conditions/results/Sign_Concordance_Summary.csv")
add("F3S_YPD1v2_sign_conc", "Supp. Fig. 3", "Sign concordance for eQTLs significant in both YPD.1 and YPD.2.", rd(f), f)
f <- file.path(F3S, "reproducibility_of_eQTL_mapping_between_conditions/results/shared_eQTLs_YPD1_vs_YPD2_ALL_variants.csv")
add("F3S_YPD1v2_shared_eQTLs", "Supp. Fig. 3", "eQTLs significant in both YPD.1 and YPD.2 (beta, -log10 P, permutation threshold).", rd(f), f)

# ---------------- Figure 4 (condition comparison) ----------------
vt4 <- c(SNPs = "snp", InDels = "indel", CNVs = "cnv", SVs = "sv")
p <- as.vector(outer(names(F4_DATASETS), vt4, function(d, t)
  file.path("Figure.4", paste0("Fig.4-", F4_DATASETS[d]), "results", F4_DATASETS[d],
            "annotations", paste0(t, "_data_annotated_cis_trans.txt"))))
ids <- as.vector(outer(names(F4_DATASETS), names(vt4), paste, sep = "|"))
add("F4_eQTL_local_distant", "Fig. 4 / Supp.",
    "All clumped eQTLs for YPD.2, Flu.1, Flu.2, sc2glu.2 and scgal.1 with local/distant (and old cis/trans) annotation. YPD.1 is sheet F3_eQTL_YPD1_cis_trans.",
    stack_tc(p, ids, "dataset_variant") %>%
      separate(dataset_variant, c("condition", "variant_type"), sep = "\\|"),
    file.path("Figure.4/Fig.4-<dataset>/results/<dataset>/annotations", "<type>_data_annotated_cis_trans.txt"))
f <- file.path(F4M, "replicates/dataset_overview_table.csv")
add("F4_dataset_overview", "Fig. 4", "Per dataset: strains, cells, UMI depth, genes tested, eQTLs, eGenes.", rd(f), f)
f <- file.path(F4M, "replicates/combined_eQTL_summary_all_datasets.csv")
add("F4_eQTL_summary_datasets", "Fig. 4", "eQTL counts, local/distant split, effect sizes, enrichment and pleiotropy per dataset and variant type.", rd(f), f)
f <- file.path(F4M, "flu_eQTLs_in_baselines/eQTLs_found_in_other_condition_summary.csv")
add("F4_detection_rates", "Fig. 4", "eQTLs of one dataset looked up in another: fraction found and effect-direction agreement (integrated datasets removed).",
    rd(f) %>% filter(!str_detect(comparison, "integrated")), f)
f <- file.path(F4M, "replicates_concordance/Concordance_summary_all_pairs.csv")
add("F4_concordance_summary", "Fig. 4", "Pairwise eQTL effect concordance between datasets (sign, correlation, permutation P).", paper_pairs(rd(f)), f)
f <- file.path(F4M, "replicates_concordance/Concordance_by_variant_type_all_pairs.csv")
add("F4_concordance_by_type", "Fig. 4", "Effect-direction concordance and correlation between datasets, by variant type.", paper_pairs(rd(f)), f)
f <- file.path(F4M, "replicates_concordance/eQTL_sharing_categories_all_pairs.csv")
add("F4_eQTL_sharing", "Fig. 4", "eQTL sharing categories (shared / only in one dataset) for each dataset pair.", paper_pairs(rd(f)), f)
f <- file.path(F4M, "replicates_concordance/Concordant_genes_summary_all_pairs.csv")
add("F4_concordant_genes_summary", "Fig. 4", "Genes with high vs low eQTL effect concordance between datasets.", paper_pairs(rd(f)), f)
f <- file.path(F4M, "replicates_concordance/Gene_concordance_stats_all_pairs.csv")
add("F4_gene_concordance", "Fig. 4", "Per-gene eQTL effect concordance for each dataset pair.", paper_pairs(rd(f)), f)
f <- file.path(F4M, "replicates_concordance/GO_gene_categories_summary_all_pairs.csv")
add("F4_gene_categories", "Fig. 4", "Genes per category (dataset-specific / shared / concordant) for each dataset pair.",
    rd(f) %>% filter(!str_detect(Comparison, "sc2glu\\.1")), f)
go_pairs <- c("YPD.1_vs_Flu.1_16ug", "YPD.1_vs_Flu.2_8ug", "YPD.1_vs_SCGal.1", "YPD.2_vs_SCGal.1",
              "SC2Glu.2_vs_SCGal.1", "SCGal.1_vs_Flu.1_16ug", "SCGal.1_vs_Flu.2_8ug")
p <- file.path(F4M, "replicates", paste0("GO_comparison_", go_pairs, "_GO_terms.csv"))
add("F4_GO_condition_specific", "Fig. 4", "GO enrichment of dataset-specific vs shared eGenes for each dataset pair.",
    stack(p, go_pairs, "comparison"), p)
f <- file.path(F4M, "mashr/dataset_qc.csv")
add("F4_mash_dataset_qc", "Fig. 4", "mashr input datasets: cells, strains, UMI depth, median posterior effect.", rd(f), f)
f <- file.path(F4M, "mashr/per_condition_summary_all_runs.csv")
add("F4_mash_per_condition", "Fig. 4", "mashr: significant and condition-specific eQTLs per dataset.", rd(f), f)
f <- file.path(F4M, "mashr/n_conditions_shared.csv")
add("F4_mash_n_shared", "Fig. 4", "mashr: number of datasets in which each eQTL is significant / shared.", rd(f), f)
p <- file.path(F4M, "mashr", c("pairwise_sharing_sign_A_main.csv", "pairwise_sharing_magnitude_A_main.csv"))
add("F4_mash_pairwise_sharing", "Fig. 4", "mashr: pairwise eQTL sharing between datasets, by sign and by magnitude.",
    stack(p, c("sign", "magnitude"), "sharing"), p)
f <- file.path(F4M, "mashr/mash_results_A_main.csv")
add("F4_mash_results", "Fig. 4", "mashr posterior effect, SD and lfsr of every eQTL in each of the 6 datasets, and condition-specific calls.", rd(f), f)
f <- file.path(F4M, "mashr/fluconazole_specific_eQTLs_robust_A_and_C.csv")
add("F4_mash_flu_specific", "Fig. 4", "Fluconazole-specific eQTLs called in both mashr runs (6 datasets; YPD + fluconazole only).", rd(f), f)
f <- file.path(F4S, "Layer1_Layer2_summary_stats.csv")
add("F4S_bulk_vs_sc_summary", "Supp. Fig. 4", "SC 2% glucose: agreement of 10x pseudobulk with published bulk RNA-seq (global and per-gene).", rd(f), f)
p <- file.path(F4S, c("Layer1b_per_gene_correlation.csv", "Layer2_per_gene_variance.csv"))
add("F4S_bulk_vs_sc_per_gene", "Supp. Fig. 4", "Per-gene bulk vs single-cell correlation and across-strain variance.",
    full_join(rd(p[1]), rd(p[2]), by = "gene"), p)

# ---------------- Figure 5 (cell cycle) ----------------
f <- file.path(F5, "harmonic_q1_pooled.csv")
add("F5_periodic_genes_pooled", "Fig. 5", "Harmonic fit of expression vs inferred cell-cycle time, all cells pooled: R2, amplitude, peak, periodicity FDR.", rd(f), f)
f <- file.path(F5, "harmonic_gene_summary.csv")
add("F5_gene_variation_summary", "Fig. 5", "Per periodic gene: across-strain variation in cell-cycle structure (R2), amplitude and peak timing vs permutation null.", rd(f), f)
f <- file.path(F5, "sig_structure.csv")
add("F5_sig_structure", "Fig. 5", "Genes whose periodic structure (R2) varies significantly between strains.", rd(f), f)
f <- file.path(F5, "sig_amplitude.csv")
add("F5_sig_amplitude", "Fig. 5", "Genes whose periodic amplitude varies significantly between strains.", rd(f), f)
f <- file.path(F5, "sig_timing.csv")
add("F5_sig_timing", "Fig. 5", "Genes whose peak timing varies significantly between strains (early/late strains listed).", rd(f), f)

# ── Write workbook ──────────────────────────────────────────────────
wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", fgFill = "#DCE6F1", border = "bottom")

index <- tibble(
  Sheet       = names(sheets),
  Figure      = map_chr(sheets, "figure"),
  Description = map_chr(sheets, "description"),
  Rows        = map_int(sheets, ~ nrow(.x$data)),
  Columns     = map_int(sheets, ~ ncol(.x$data)),
  Source      = map_chr(sheets, "sources")
)

addWorksheet(wb, "README")
writeData(wb, "README", "Supplementary tables — Agrawal et al. single-cell eQTL / cell cycle. One sheet per table; click a sheet name to jump to it.",
          startRow = 1)
writeData(wb, "README", paste("Source paths are relative to", PROJ), startRow = 2)
writeData(wb, "README", index, startRow = 4, headerStyle = hdr)
for (i in seq_len(nrow(index))) {
  writeFormula(wb, "README", startRow = 4 + i, startCol = 1,
               x = makeHyperlinkString(sheet = index$Sheet[i], row = 1, col = 1, text = index$Sheet[i]))
}
setColWidths(wb, "README", cols = 1:6, widths = c(30, 20, 90, 8, 8, 100))
freezePane(wb, "README", firstActiveRow = 5)

for (nm in names(sheets)) {
  d <- sheets[[nm]]$data
  addWorksheet(wb, nm)
  writeData(wb, nm, d, headerStyle = hdr, withFilter = TRUE)
  freezePane(wb, nm, firstRow = TRUE)
  setColWidths(wb, nm, cols = seq_len(ncol(d)), widths = pmin(pmax(nchar(names(d)) + 2, 10), 40))
  cat(sprintf("%-30s %7d rows\n", nm, nrow(d)))
}

saveWorkbook(wb, OUT, overwrite = TRUE)
cat("\nWrote", OUT, "with", length(sheets), "data sheets + README\n")
