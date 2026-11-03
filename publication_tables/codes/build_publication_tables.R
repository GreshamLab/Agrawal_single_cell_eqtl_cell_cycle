# ============================================================
# Publication supplementary tables, numbered in the order they are cited in
# the Results text (Supplementary Tables 1-22, plus main Table 1).
#
# Output (publication_tables/results/):
#   Agrawal_Supplementary_Tables.xlsx  README index + one sheet per table
#   csv/<sheet>.csv                    the same tables as CSV
#
# Every table is read from the figure folders (paths in the README sheet);
# nothing is recomputed except where noted (S1 CNV subtypes, S13B random set).
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
})

PROJ <- "/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle"
OUT_DIR <- file.path(PROJ, "publication_tables/results")
dir.create(file.path(OUT_DIR, "csv"), recursive = TRUE, showWarnings = FALSE)

F1   <- "Figure.1"
F2S  <- "Figure.2/Fig.2_supplementary/results"
F3M  <- "Figure.3/Fig.3-Panel-A-C-D-cis-and-trans-eqtl-YPD-variant-type/results"
F4M  <- "Figure.4/Fig.4-PartII-compare-conditions/results"
F5   <- "Figure.5/2.cell_cycle_time_assignment_to_single_cells/results"
F5T  <- file.path(F5, "paper_tables")

# The six paper datasets (sc2glu.1 and the integrated datasets are not in the tables).
# Flu.1 = 16 ug/ml and Flu.2 = 8 ug/ml: two doses, not replicates.
DATASETS <- tribble(
  ~key,       ~label,                                ~dir,
  "YPD.1",    "YPD.1",                               file.path(F3M, "YPD.1"),
  "YPD.2",    "YPD.2",                               "Figure.4/Fig.4-YPD.2/results/YPD.2",
  "Flu.1",    "Fluconazole 16 ug/ml (Flu.1)",        "Figure.4/Fig.4-Flu.1/results/Flu.1",
  "Flu.2",    "Fluconazole 8 ug/ml (Flu.2)",         "Figure.4/Fig.4-Flu.2/results/Flu.2",
  "sc2glu.2", "SC glucose (sc2glu.2)",               "Figure.4/Fig.4-sc2glu.2_only/results/sc2glu.2_only",
  "scgal.1",  "SC galactose (scgal.1)",              "Figure.4/Fig.4-scgal.1/results/scgal.1"
)
ds_label <- setNames(DATASETS$label, DATASETS$key)

rd <- function(path, ...) {
  full <- file.path(PROJ, path)
  if (grepl("\\.csv$", path)) read_csv(full, show_col_types = FALSE, guess_max = 1e5, ...)
  else read_tsv(full, show_col_types = FALSE, guess_max = 1e5, ...)
}
# Read one file per dataset and stack them with a leading Dataset column
stack_ds <- function(rel) {
  map2(DATASETS$dir, DATASETS$label, ~ rd(file.path(.x, rel)) %>% mutate(Dataset = .y, .before = 1)) %>%
    bind_rows()
}

# ── Sheet registry ─────────────────────────────────────────────────
sheets <- list()
add <- function(sheet, table, section, description, data, sources) {
  stopifnot(nchar(sheet) <= 31, !sheet %in% names(sheets))
  sheets[[sheet]] <<- list(table = table, section = section, description = description,
                           data = data, sources = paste(sources, collapse = "; "))
}

# ================= Figure 1 =================
# S1. Variants tested (MAF > 15%). SNP/InDel/SV rows from the Fig 1 table; gene-based
# CNV rows recomputed from the .bim used for eQTL mapping, with the same filters as
# the eQTL pipeline (chr17 removed; CNVs past the S288C chromosome end removed),
# so the total matches the 3,062 CNVs tested.
f_t1 <- file.path(F1, "Fig.1-all-experiment-table/results/Table1_variant_type_subtype_summary.csv")
cnv_bim <- "/projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/structural_variant_vcf/SVcalls2/SVcalls/CNVs_96_strains/CNV_96strains_maf15/no_chr17/CNVs.96strains.maf15.no_chr17.plink.bim"
chr_lengths <- c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643,
                 439888, 745751, 666816, 1078177, 924431, 784333, 1091291, 948066)
cnv <- read.table(cnv_bim, header = FALSE, stringsAsFactors = FALSE) %>%
  filter(V1 != 17, V4 <= chr_lengths[V1]) %>%
  mutate(Subtype = ifelse(as.numeric(sub("^.*_Cov", "", V2)) < 1, "Loss", "Gain"))
s1 <- rd(f_t1) %>%
  filter(`Variant Type` != "Gene-based CNVs") %>%
  bind_rows(count(cnv, Subtype, name = "Total number of variants") %>%
              mutate(`Variant Type` = "Gene-based CNVs", .before = 1))
s1 <- bind_rows(s1, s1 %>% group_by(`Variant Type`) %>%
                  summarise(Subtype = "All", `Total number of variants` = sum(`Total number of variants`))) %>%
  arrange(factor(`Variant Type`, c("SNPs", "InDels", "SVs", "Gene-based CNVs")), Subtype == "All")
add("S1_variants", "Supp. Table 1", "Fig. 1",
    "Common variants (MAF > 15%) in the 96-strain panel by type and subtype (SNPs 51,796; InDels 3,021; SVs 258; gene-based CNVs 3,062).",
    s1, c(f_t1, paste(cnv_bim, "(CNV rows; chr17 and positions past chromosome end removed)")))

f <- file.path(F1, "Fig.1-all-experiment-table/results/Table2_pairwise_diversity_summary.csv")
add("S2_pairwise_diversity", "Supp. Table 2", "Fig. 1",
    "Pairwise differences between the 96 strains per variant type (median 25,314 SNPs, 1,463 InDels, 116 SVs, 1,243 gene-based CNVs).",
    rd(f), f)

# S3. Cells per run, with demuxlet singlet/doublet/ambiguous calls
f <- file.path(F1, "Fig.1-all-experiment-table/results/Supplementary_Table_QC_summary.csv")
run_label <- c("YPD.1" = "YPD.1", "YPD.2" = "YPD.2",
               "Fluconazole.1" = "Fluconazole 16 ug/ml (Flu.1)", "Fluconazole.2" = "Fluconazole 8 ug/ml (Flu.2)",
               "Synthetic-complete.1" = "SC glucose replicate 1 (sc2glu.1)",
               "Synthetic-complete.2" = "SC glucose replicate 2 (sc2glu.2)",
               "Galactose.1" = "SC galactose (scgal.1)")
s3 <- rd(f) %>%
  mutate(Dataset = unname(run_label[paste(Condition, Replicate, sep = ".")]), .before = 1,
         Total_reads = as.numeric(gsub(",", "", Total_reads)))
s3 <- bind_rows(s3, s3 %>% summarise(Dataset = "Total", across(c(N_cells_total, Total_reads,
                                                                  N_singlets_demuxlet, N_doublets_demuxlet,
                                                                  N_ambiguous_demuxlet), sum)))
add("S3_cells_per_dataset", "Supp. Table 3", "Fig. 1",
    "Cells, reads and genes per scRNA-seq run, demuxlet singlet/doublet/ambiguous calls (257,249 cells in total), and the number of strains demuxlet assigned singlets to before (N_strains_demuxlet) and after (N_strains_min5_cells) the >= 5 cells per strain filter.",
    s3, f)

# S4. Nested PERMANOVA (condition vs replicate run within condition). adonis2 output is
# only printed in the knitted HTML, so it is parsed from there.
f_html <- file.path(F1, "Fig.1-Panel-C-scRNAseq-reproducibility_agreement/codes/Fig1-PanelC-does-gene-experession-differ-bwn-condition.html")
html <- readLines(file.path(PROJ, f_html), warn = FALSE)
parse_adonis <- function(obj, analysis) {
  i <- grep(paste0("adonis2(formula = dist(", obj, "$x) ~ condition/replicate"), html, fixed = TRUE)[1]
  stopifnot(!is.na(i))
  rows <- html[(i + 2):(i + 5)] %>% sub("^## *", "", .) %>% sub("</code>.*", "", .) %>% trimws()
  map_dfr(rows, function(r) {
    term <- str_extract(r, "^\\S+")
    v <- str_split(trimws(sub("^\\S+", "", r)), "\\s+")[[1]]
    tibble(Term = term, Df = as.numeric(v[1]), SumOfSqs = as.numeric(v[2]), R2 = as.numeric(v[3]),
           F = suppressWarnings(as.numeric(v[4])), P = suppressWarnings(as.numeric(v[5])))
  }) %>% mutate(Analysis = analysis, .before = 1)
}
s4a <- bind_rows(
  parse_adonis("res_all_gal", "All datasets + SC galactose (7 datasets, 58 strains, 500 genes)"),
  parse_adonis("res_all", "All datasets without galactose (6 datasets, 64 strains, 500 genes)")
) %>% mutate(Term = recode(Term, "condition:replicate" = "replicate run within condition"),
             Pct_variance = round(100 * R2, 2), .after = R2)
add("S4A_PERMANOVA", "Supp. Table 4", "Fig. 1E",
    "Nested PERMANOVA (adonis2, Euclidean, 999 permutations) of strain pseudobulk profiles: condition 51.3% vs replicate run within condition 2.2% of variance.",
    s4a, paste(f_html, "(printed adonis2 output)"))
p <- file.path(F1, "Fig.1-Panel-C-scRNAseq-reproducibility_agreement/results",
               c("Fig1C_condition_PCA_with_galactose_summary_stats.csv", "Fig1C_condition_PCA_summary_stats.csv"))
add("S4B_PCA_summary", "Supp. Table 4", "Fig. 1E",
    "PCA of strain pseudobulk profiles per dataset set: PC1/PC2 variance, one-factor PERMANOVA R2 for condition, silhouette.",
    map_dfr(p, rd), p)

# ================= Figure 2 =================
p <- file.path(F2S, DATASETS$key, "ora_tables", paste0("clade_summary_", DATASETS$key, ".csv"))
add("S5A_DEGs_per_clade", "Supp. Table 5", "Fig. 2",
    "Differentially expressed genes per clade (clade vs rest, DESeq2 on pseudobulk; padj < 0.05, |log2FC| > 0.3), each dataset.",
    map_dfr(p, rd) %>% mutate(condition = recode(condition, !!!ds_label)), p)

p_de <- file.path(F2S, DATASETS$key, "ora_tables", paste0("DESeq2_clade_vs_rest_all_genes_", DATASETS$key, ".csv"))
de <- setNames(map(p_de, rd), DATASETS$key)
add("S5B_unique_DEGs_per_dataset", "Supp. Table 5", "Fig. 2",
    "Unique differentially expressed genes per dataset (padj < 0.05, |log2FC| > 0.3 in at least one clade-vs-rest contrast).",
    imap_dfr(de, ~ .x %>% filter(padj < 0.05, abs(log2FoldChange) > 0.3) %>%
               summarise(Dataset = ds_label[.y], clades_tested = n_distinct(.x$clade),
                         unique_DEGs = n_distinct(gene),
                         unique_upregulated = n_distinct(gene[log2FoldChange > 0]),
                         unique_downregulated = n_distinct(gene[log2FoldChange < 0]))),
    p_de)
p <- file.path(F2S, DATASETS$key, "ora_tables", paste0("GO_ORA_", DATASETS$key, ".csv"))
add("S5C_GO_ORA_per_clade", "Supp. Table 5", "Fig. 2B-D",
    "GO over-representation (clusterProfiler) of up- and down-regulated DEGs of each clade, each dataset.",
    map2_dfr(p, DATASETS$label, ~ rd(.x) %>% mutate(Dataset = .y, .before = 1)), p)
for (i in seq_len(nrow(DATASETS))) {
  add(paste0("S5", LETTERS[3 + i], "_DESeq2_", DATASETS$key[i]), "Supp. Table 5", "Fig. 2",
      paste0("Full DESeq2 clade-vs-rest results for every tested gene, ", DATASETS$label[i],
             " (design ~group; DEG = padj < 0.05 & |log2FC| > 0.3)."),
      de[[i]], p_de[i])
}

# ================= Figure 3 / 4: per-dataset eQTLs =================
rel_eqtl <- "eqtl_plot/annotations/cis_trans_eQTL_mapping_all_variants_combined.txt"
for (i in seq_len(nrow(DATASETS))) {
  f <- file.path(DATASETS$dir[i], rel_eqtl)
  add(paste0("S6", LETTERS[i], "_eQTLs_", DATASETS$key[i]), "Supp. Table 6", "Fig. 3 / Fig. 4",
      paste0("All LD-clumped eQTLs (SNPs, InDels, SVs, gene-based CNVs) in ", DATASETS$label[i],
             ": target gene, variant, causal gene, local/distant, |beta|, P."),
      rd(f), f)
}

# S7. Single-cell NB GLM re-test of YPD.1 local eQTLs
vt <- c("SNPs", "InDels", "SVs", "CNVs")
p <- file.path(F3M, "raw_single_cell_plots", paste0("PanelB_verify_", vt, "_local.csv"))
add("S7A_local_eQTL_singlecell", "Supp. Table 7", "Fig. 3C",
    "YPD.1 local eQTLs re-tested in single cells: negative binomial GLM of per-cell UMI counts on strain genotype vs FaST-LMM beta/P.",
    map2_dfr(p, vt, ~ rd(.x) %>% mutate(Variant_type = .y, .before = 1)), p)
f <- file.path(F3M, "raw_single_cell_plots/PanelB_pvalue_and_effect_agreement_stats.csv")
add("S7B_local_eQTL_agreement", "Supp. Table 7", "Fig. 3C",
    "Summary: direction concordance (exact binomial test) and Spearman correlation of single-cell vs FaST-LMM effects.",
    rd(f), f)

# Main Table 1 and S8-S12: per-dataset eQTL summaries
f <- file.path(F4M, "replicates/dataset_overview_table.csv")
add("Table1_dataset_overview", "Main Table 1", "Fig. 4",
    "Per dataset: strains, cells, UMI depth, genes tested, between-strain expression variance, eQTLs, eGenes.",
    rd(f), f)
f_sum <- file.path(F4M, "replicates/combined_eQTL_summary_all_datasets.csv")
es <- rd(f_sum) %>% mutate(Dataset = recode(key, !!!ds_label), .before = 1) %>% select(-dataset, -key)
add("S8_eQTLs_per_dataset", "Supp. Table 8", "Fig. 4",
    "eQTLs per dataset after per-gene LD clumping: associations, unique eQTL variants, eGenes.",
    es %>% group_by(Dataset, n_strains, n_genes_tested) %>%
      summarise(eQTL_associations = sum(eQTLs), unique_eQTL_variants = sum(eQTL_variants),
                local = sum(local), distant = sum(distant), .groups = "drop") %>%
      left_join(rd(f) %>% transmute(key = c("YPD.1", "YPD.2", "sc2glu.2", "Flu.1", "Flu.2", "scgal.1")[match(dataset, c("YPD.1", "YPD.2", "synthetic complete", "fluconazole 16 ug/ml", "fluconazole 8 ug/ml", "galactose"))],
                                    Dataset = ds_label[key], eGenes) %>% select(-key), by = "Dataset") %>%
      arrange(match(Dataset, DATASETS$label)),
    c(f_sum, f))
add("S9_eQTL_variant_types", "Supp. Table 9", "Fig. 4",
    "Unique eQTL variants by variant type per dataset, as a share of eQTL variants and of variants tested.",
    es %>% group_by(Dataset) %>%
      transmute(Genetic_variants, tested_variants, eQTL_variants,
                pct_of_eQTL_variants = round(100 * eQTL_variants / sum(eQTL_variants), 1),
                pct_of_tested_that_are_eQTL = round(100 * eQTL_variants / tested_variants, 1)) %>% ungroup(),
    f_sum)
add("S10_local_distant", "Supp. Table 10", "Fig. 4",
    "Local vs distant eQTLs per dataset and variant type, with median |beta| of each class.",
    es %>% select(Dataset, Genetic_variants, eQTLs, local, distant, pct_distant,
                  median_abs_beta, median_abs_beta_local, median_abs_beta_distant),
    f_sum)
rel <- "eqtl_plot/data_summary/odds_ratio_enrichment.txt"
add("S11_eQTL_enrichment", "Supp. Table 11", "Fig. 3F / Fig. 4",
    "Fisher's exact test of each variant type among eQTL variants vs all other types (a = eQTL, b = not eQTL; c/d = other types), FDR-adjusted.",
    stack_ds(rel), file.path(DATASETS$dir, rel))
rel <- "eqtl_plot/data_summary/pleiotropy_summary.txt"
add("S12_pleiotropy", "Supp. Table 12", "Fig. 4",
    "Pleiotropy: genes affected per eQTL variant, by variant type and dataset.",
    stack_ds(rel), file.path(DATASETS$dir, rel))

# ================= mashr =================
MASH <- file.path(F4M, "mashr")
strong <- rd(file.path(MASH, "mash_results_A_main.csv"))
# Random set: the cached random pairs that mashr used, i.e. those with a finite
# effect and SE > 0 in all six datasets (as in make_mats() of the mashr Rmd)
ex_names <- c("YPD.1", "YPD.2", "Flu_16ug", "Flu_8ug", "sc2glu.2", "galactose")
random <- readRDS(file.path(PROJ, MASH, "random_pairs.rds"))
for (cn in ex_names) {
  e <- readRDS(file.path(PROJ, MASH, paste0("extract_", cn, ".rds"))) %>%
    distinct(Gene, SNP, Variant_Type, .keep_all = TRUE) %>%
    select(Gene, SNP, Variant_Type, beta = SnpWeight, se = SnpWeightSE)
  random <- left_join(random, e, by = c("Gene", "SNP", "Variant_Type")) %>%
    rename(!!paste0("beta_", cn) := beta, !!paste0("se_", cn) := se)
}
random <- random %>%
  filter(if_all(starts_with("beta_"), is.finite), if_all(starts_with("se_"), ~ is.finite(.x) & .x > 0))
add("S13A_mash_set_sizes", "Supp. Table 13", "Fig. 4",
    "mashr input sets by variant type: strong set (lead variant per gene) and random set (random tested pairs).",
    full_join(count(strong, Variant_Type, name = "strong_set"), count(random, Variant_Type, name = "random_set"),
              by = "Variant_Type") %>%
      bind_rows(summarise(., Variant_Type = "All", across(-Variant_Type, sum))),
    c(file.path(MASH, "mash_results_A_main.csv"), file.path(MASH, "random_pairs.rds")))
add("S13B_mash_strong_set", "Supp. Table 13", "Fig. 4",
    "Strong set (4,973 lead eQTLs, one per gene; eQTLs without a local / distant call excluded): mashr posterior effect, SD and lfsr in each dataset, number of datasets significant, condition-specific calls.",
    strong, file.path(MASH, "mash_results_A_main.csv"))
add("S13C_mash_random_set", "Supp. Table 13", "Fig. 4",
    "Random set (16,885 gene-variant pairs complete in all six datasets; reproducible draw, seed 42): FaST-LMM effect and SE in each dataset.",
    random, c(file.path(MASH, "random_pairs.rds"), file.path(MASH, "extract_<dataset>.rds")))

# ================= Figure 5 =================
genes <- rd(file.path(F5T, "SuppTable_gene_phenotypes.csv"))
f <- file.path(F5T, "SuppTable_gene_phenotypes.csv")
add("S14_cell_cycle_gene_variation", "Supp. Table 14", "Fig. 5",
    "All 478 Spellman cell-cycle genes: per-strain harmonic fit summary, between-strain amplitude and timing spread vs permutation null, FDR and significance.",
    genes, f)
f <- file.path(F5T, "Table5B_summary.csv")
add("S15_variable_gene_counts", "Supp. Table 15", "Fig. 5",
    "Genes with significant between-strain variation in amplitude (168), timing (101), either (211) or both (58).",
    rd(f), f)
f <- file.path(F5T, "SuppTable_histone_robustness.csv")
add("S16_histone_robustness", "Supp. Table 16", "Fig. 5C",
    "MRH1, TIP1 and histone genes: spread after trimming the 10% most extreme strains, robustness, and effect of excluding strain APQ.",
    rd(f), f)
f <- file.path(F5T, "SuppTable_gene_phenotypes.csv")
add("S17_amplitude_ranking", "Supp. Table 17", "Fig. 5",
    "478 genes ranked by between-strain amplitude variation (trimmed spread relative to permutation null).",
    genes %>% arrange(desc(`Amplitude trimmed excess`)) %>%
      select(Gene, ORF, `Amplitude trimmed excess`, `Amplitude robustness`, `Amplitude CV across strains`,
             `Amplitude CV, permutation null`, `Fraction of amplitude spread explained by expression level`,
             `Amplitude FDR`, `Significant amplitude variation`),
    f)
add("S18_timing_ranking", "Supp. Table 18", "Fig. 5",
    "478 genes ranked by between-strain peak-time variation (trimmed spread relative to permutation null).",
    genes %>% arrange(desc(`Timing trimmed excess`)) %>%
      select(Gene, ORF, `Timing trimmed excess`, `Timing robustness`, `Mean peak time (min)`,
             `Peak time circular SD (min)`, `Peak time SD, permutation null (min)`, `Timing FDR`,
             `Significant timing variation`, `Earliest strains (up to 5)`, `Latest strains (up to 5)`),
    f)
f_go <- file.path(F5T, "SuppTable_GO_enrichment_full.csv")
go <- rd(f_go) %>% filter(background == "Spellman genes") %>% arrange(set, p)
add("S19_GO_variable_genes", "Supp. Table 19", "Fig. 5",
    "GO enrichment of amplitude-variable (168) and timing-variable (101) genes in YPD.1, background = 478 Spellman genes.",
    go %>% filter(set %in% c("Amplitude-variable (YPD.1)", "Timing-variable (YPD.1)")), f_go)
f <- file.path(F5T, "SuppTable_replicate_overlap.csv")
add("S20_replicate_overlap", "Supp. Table 20", "Fig. 5 replicate",
    "Genes significant for amplitude or timing variation in YPD.1, YPD.2 and both (47 amplitude, 21 timing reproduced).",
    rd(f), f)
f <- file.path(F5T, "SuppTable_replicate_reproducibility.csv")
add("S21_replicate_rank_corr", "Supp. Table 21", "Fig. 5 replicate",
    "Per gene: significance in each replicate and Spearman correlation (rho, FDR) of per-strain amplitude and peak time between YPD.1 and YPD.2.",
    rd(f), f)
add("S22_GO_reproduced_genes", "Supp. Table 22", "Fig. 5 replicate",
    "GO enrichment of genes whose amplitude (47) or timing (21) variation reproduced in YPD.2, background = 478 Spellman genes.",
    go %>% filter(set %in% c("Amplitude reproduced in YPD.2", "Timing reproduced in YPD.2")), f_go)

# ── Write workbook + CSVs ──────────────────────────────────────────
wb  <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", fgFill = "#DCE6F1", border = "bottom")

index <- tibble(
  Table       = map_chr(sheets, "table"),
  Sheet       = names(sheets),
  Figure      = map_chr(sheets, "section"),
  Description = map_chr(sheets, "description"),
  Rows        = map_int(sheets, ~ nrow(.x$data)),
  Source      = map_chr(sheets, "sources")
)

addWorksheet(wb, "README")
writeData(wb, "README", "Agrawal et al. single-cell eQTL / cell cycle: supplementary tables, numbered in the order cited in the Results. Click a sheet name to open it.", startRow = 1)
writeData(wb, "README", paste("Source paths are relative to", PROJ, "unless absolute."), startRow = 2)
writeData(wb, "README", index, startRow = 4, headerStyle = hdr)
for (i in seq_len(nrow(index))) {
  writeFormula(wb, "README", startRow = 4 + i, startCol = 2,
               x = makeHyperlinkString(sheet = index$Sheet[i], row = 1, col = 1, text = index$Sheet[i]))
}
setColWidths(wb, "README", cols = 1:6, widths = c(15, 30, 16, 90, 8, 110))
freezePane(wb, "README", firstActiveRow = 5)

for (nm in names(sheets)) {
  d <- sheets[[nm]]$data
  addWorksheet(wb, nm)
  writeData(wb, nm, d, headerStyle = hdr, withFilter = TRUE)
  freezePane(wb, nm, firstRow = TRUE)
  setColWidths(wb, nm, cols = seq_len(ncol(d)), widths = pmin(pmax(nchar(names(d)) + 2, 10), 40))
  write_csv(d, file.path(OUT_DIR, "csv", paste0(nm, ".csv")), na = "")
  cat(sprintf("%-32s %7d rows\n", nm, nrow(d)))
}
write_csv(index, file.path(OUT_DIR, "README_table_index.csv"))

out <- file.path(OUT_DIR, "Agrawal_Supplementary_Tables.xlsx")
saveWorkbook(wb, out, overwrite = TRUE)
cat("\nWrote", out, "with", length(sheets), "sheets + README\n")
