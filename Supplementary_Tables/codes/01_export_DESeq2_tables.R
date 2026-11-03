# ============================================================
# Export full per-clade DESeq2 results (all genes) for every condition.
#
# The Figure 2 Rmds compute these (dge_combined) but only save GO/GSEA
# and DEG counts. This script reproduces the identical DE step —
# same pseudobulk inputs, ORF conversion, chromosome filter, clade
# metadata, and run_clade_dge() (focal clade vs all other strains,
# design ~ group, padj < 0.05 & |log2FC| > 0.3) — and writes the tables.
#
# Sources:
#   Figure.2/Fig.2-Panel-A-B-C-D-.../codes/Fig.2-panel-A-B-C-D-differential-gene-expression.Rmd
#   Figure.2/Fig.2_supplementary/codes/fig.2.supplementary-replicates-per-condition-DEG.Rmd
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(org.Sc.sgd.db)
  library(AnnotationDbi)
  library(parallel)
})

select <- dplyr::select
filter <- dplyr::filter
first  <- dplyr::first

BASE    <- "/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all"
PROJ    <- file.path(BASE, "Agrawal_single_cell_eqtl_cell_cycle")
ANNOT   <- "/scratch/aa9679/SingleCell_cerevisiae/Gene_exp_analysis/yeast_custom_annotations.csv"
OUT_DE  <- file.path(PROJ, "Supplementary_Tables/results/DESeq2_tables")
FIG2_MAIN_OUT <- file.path(PROJ, "Figure.2/Fig.2-Panel-A-B-C-D-strain-clade-gene-expression-comparison/results")
FIG2_SUPP_OUT <- file.path(PROJ, "Figure.2/Fig.2_supplementary/results")
dir.create(OUT_DE, showWarnings = FALSE, recursive = TRUE)

N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

CONDITIONS <- list(
  "YPD.1"    = list(path = file.path(BASE, "YPD.1/fastlmm/data/pheno_data/pseudobulk_reliable_transcripts_YPD.txt"),    sep = "\t"),
  "YPD.2"    = list(path = file.path(BASE, "YPD.2/fastlmm/data/pheno_data/pseudobulk_reliable_transcripts_YPD.txt"),    sep = "\t"),
  "Flu.1"    = list(path = file.path(BASE, "Flu.1/fastlmm/data/pheno_data/pseudobulk_reliable_transcripts_Flu1.txt"),   sep = "\t"),
  "Flu.2"    = list(path = file.path(BASE, "Flu.2/fastlmm/data/pheno_data/pseudobulk_reliable_transcripts_Flu2.txt"),   sep = "\t"),
  "sc2glu.2" = list(path = file.path(BASE, "sc2glu.2/fastlmm/data/pheno_data/pseudobulked_sc2glu_before_any_filtering.csv"), sep = ","),
  "scgal.1"  = list(path = file.path(BASE, "scgal.1/fastlmm/data/pheno_data/pseudobulk_reliable_transcripts_scgal.1.txt"), sep = "\t")
)

# ── Shared annotation resources ─────────────────────────────────────
all_custom_gene_annotations <- read.csv(ANNOT)

common_to_orf <- all_custom_gene_annotations %>%
  filter(!is.na(gene_name), !is.na(gene_id)) %>%
  distinct(gene_name, gene_id) %>%
  deframe()

orf_to_name_custom <- all_custom_gene_annotations %>%
  filter(!is.na(gene_id), !is.na(gene_name)) %>%
  distinct(gene_id, gene_name) %>%
  deframe()

gene_location <- all_custom_gene_annotations %>%
  filter(type == "CDS", !is.na(gene_id)) %>%
  group_by(gene_id) %>%
  summarise(gene_name = first(gene_name), chr = first(seqnames), .groups = "drop")

standard_orfs <- gene_location %>% filter(!chr %in% c("Mito", "Marker")) %>% pull(gene_id)

missing_orfs <- unique(c(names(orf_to_name_custom), standard_orfs))
missing_orfs <- missing_orfs[!missing_orfs %in% names(orf_to_name_custom)]
sgd <- tryCatch(
  AnnotationDbi::select(org.Sc.sgd.db, keys = missing_orfs, columns = "GENENAME", keytype = "ORF"),
  error = function(e) data.frame(ORF = character(), GENENAME = character())
)
sgd_map <- sgd %>% filter(!is.na(GENENAME)) %>% distinct(ORF, GENENAME) %>% deframe()
orf_to_name <- c(orf_to_name_custom, sgd_map[!names(sgd_map) %in% names(orf_to_name_custom)])

strains_sc <- read.csv(file.path(BASE, "StrainsforSC-RNAseq.csv"), header = TRUE)

# ── Helpers (identical logic to the Figure 2 Rmds) ──────────────────
load_condition_counts <- function(path, sep = "\t") {
  quote_arg <- if (sep == ",") "\"" else ""
  raw_counts <- read.table(path, header = TRUE, sep = sep, row.names = 1,
                           check.names = FALSE, quote = quote_arg)

  expr_gene_ids <- rownames(raw_counts)
  is_orf <- grepl("^Y[A-Z]{2}\\d{3}[WC](-[A-Z])?$", expr_gene_ids)
  expr_orf_ids <- ifelse(
    is_orf, expr_gene_ids,
    ifelse(expr_gene_ids %in% names(common_to_orf), common_to_orf[expr_gene_ids], expr_gene_ids)
  )
  not_converted <- expr_gene_ids[!is_orf & !expr_gene_ids %in% names(common_to_orf)]
  rownames(raw_counts) <- expr_orf_ids
  raw_counts <- raw_counts[!rownames(raw_counts) %in% not_converted, ]
  raw_counts <- raw_counts[rownames(raw_counts) %in% standard_orfs, ]
  round(as.matrix(raw_counts))
}

build_sample_meta <- function(counts_mat) {
  strains_sc %>%
    select(strain = Standardized.name, clade = Clades) %>%
    mutate(clade = str_trim(clade)) %>%
    filter(strain %in% colnames(counts_mat)) %>%
    distinct(strain, clade) %>%
    arrange(match(strain, colnames(counts_mat))) %>%
    column_to_rownames("strain")
}

run_clade_dge <- function(clade_of_interest, counts_mat, meta,
                          padj_thresh = 0.05, lfc_thresh = 0.3) {
  n_focal <- sum(meta$clade == clade_of_interest, na.rm = TRUE)
  if (n_focal < 2) return(NULL)

  meta_de <- meta %>%
    mutate(group = factor(ifelse(clade == clade_of_interest, "focal", "rest"),
                          levels = c("rest", "focal")))
  dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = meta_de, design = ~ group)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", "focal", "rest"), alpha = padj_thresh)

  as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    mutate(
      gene_name = ifelse(gene %in% names(orf_to_name), orf_to_name[gene], gene),
      direction = case_when(
        padj < padj_thresh & log2FoldChange >  lfc_thresh ~ "Upregulated",
        padj < padj_thresh & log2FoldChange < -lfc_thresh ~ "Downregulated",
        TRUE ~ "NS"),
      clade = clade_of_interest, n_focal = n_focal)
}

# ── Run all conditions ──────────────────────────────────────────────
all_conditions <- imap(CONDITIONS, function(cfg, label) {
  cat("\n==== ", label, " ====\n", sep = "")
  counts <- load_condition_counts(cfg$path, cfg$sep)
  meta   <- build_sample_meta(counts)
  counts <- counts[, rownames(meta), drop = FALSE]
  stopifnot(identical(rownames(meta), colnames(counts)))

  clades <- sort(unique(na.omit(meta$clade)))
  res <- mclapply(clades, run_clade_dge, counts_mat = counts, meta = meta,
                  mc.cores = N_CORES)
  dge <- bind_rows(Filter(Negate(is.null), res)) %>%
    transmute(condition = label,
              clade, n_focal_strains = n_focal,
              n_rest_strains = ncol(counts) - n_focal,
              gene, gene_name,
              baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, direction) %>%
    arrange(clade, padj)

  cat(ncol(counts), "strains |", n_distinct(dge$clade), "clades tested |",
      sum(dge$direction != "NS"), "significant clade-gene pairs |",
      n_distinct(dge$gene[dge$direction != "NS"]), "unique DEGs\n")

  write_csv(dge, file.path(OUT_DE, paste0("DESeq2_clade_vs_rest_all_genes_", label, ".csv")))
  # Also drop a copy next to the existing per-condition Figure 2 tables
  write_csv(dge, file.path(FIG2_SUPP_OUT, label, "ora_tables",
                           paste0("DESeq2_clade_vs_rest_all_genes_", label, ".csv")))
  if (label == "YPD.1")
    write_csv(dge, file.path(FIG2_MAIN_OUT, "Fig2_clade_DESeq2_all_genes_YPD1.csv"))
  dge
})

cat("\nDone.\n")
sessionInfo()
