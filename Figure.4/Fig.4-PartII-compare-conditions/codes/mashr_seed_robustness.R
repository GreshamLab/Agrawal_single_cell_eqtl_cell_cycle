# ============================================================
# Robustness of the main mashr analysis (run A, six datasets) to the random set.
# Fits mashr exactly as Fig4.mashr-condition-specific-eQTLs.Rmd (sections 2-4) with the random
# pairs drawn with one seed, and saves the fit and per-eQTL calls to
# results/mashr/seed_robustness/seed_<seed>.rds. Seed 42 is the seed of the main analysis, so
# its fit must reproduce the main fit (checked in mashr_seed_robustness_summary.R).
# Usage: Rscript mashr_seed_robustness.R <seed>
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(parallel); library(mashr)
})
select <- dplyr::select; filter <- dplyr::filter

SEED <- as.integer(commandArgs(TRUE)[1])
ROOT <- "/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all"
MASH_DIR <- file.path(ROOT, "Agrawal_single_cell_eqtl_cell_cycle/Figure.4/Fig.4-PartII-compare-conditions/results/mashr")
OUT_DIR <- file.path(MASH_DIR, "seed_robustness"); dir.create(OUT_DIR, showWarnings = FALSE)
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))

# ---- datasets: same as the main analysis ----
SNP_DIR   <- "results_SNPs_96_MAF15_from_969_strains_with_double_id_final_filtered_for_ncount_less_than_5"
INDEL_DIR <- "results_Indels_96_MAF15_with_double_id_final_filtered_for_ncount_less_than_5"
CNV_DIR   <- "results_CNVs_with_double_id_is_not_needed_final_filtered_for_ncount_less_than_5_run_from_new_plink_files_no_chr17"
SV_DIR    <- "results_SVs_with_double_id_is_not_needed_final_filtered_for_ncount_less_than_5_run_from_new_plink_files"
VARIANTS  <- c("SNPs", "InDels", "SVs", "CNVs")
make_cond <- function(res, tag, indel_dir = INDEL_DIR,
                      indel_file = paste0("all_signif_Indels_", tag, "_96_MAF15_with_double_id_final_filtered_for_ncount_less_than_5_with_threshold.txt")) {
  res <- file.path(ROOT, res)
  list(signif = c(
         SNPs   = file.path(res, paste0("all_signif_SNPs_", tag, "_96_MAF15_from_969_strains_with_double_id_final_filtered_for_ncount_less_than_5_with_threshold.txt")),
         InDels = file.path(res, indel_file),
         SVs    = file.path(res, paste0("all_signif_SVs_", tag, "_with_double_id_is_not_needed_final_filtered_for_ncount_less_than_5_run_from_new_plink_files_with_threshold.txt")),
         CNVs   = file.path(res, paste0("all_signif_CNV_", tag, "_with_double_id_is_not_needed_final_filtered_for_ncount_less_than_5_run_from_new_plink_files_no_chr17_with_threshold.txt"))),
       dirs = c(SNPs = file.path(res, SNP_DIR), InDels = file.path(res, indel_dir),
                SVs = file.path(res, SV_DIR), CNVs = file.path(res, CNV_DIR)))
}
CONDS <- list(
  YPD.1     = make_cond("YPD.1/fastlmm/results", "YPD.1",
                        indel_dir  = "results_Indels_with_double_id_final_filtered_for_ncount_less_than_5",
                        indel_file = "all_signif_Indels_YPD.1_with_double_id_final_filtered_for_ncount_less_than_5_with_threshold.txt"),
  YPD.2     = make_cond("YPD.2/fastlmm/results", "YPD.2"),
  Flu_16ug  = make_cond("Flu.1/fastlmm/results", "Flu.1"),
  Flu_8ug   = make_cond("Flu.2/fastlmm/results", "Flu.2"),
  sc2glu.2  = make_cond("sc2glu.2/fastlmm_sc2.glu.2_only/results", "sc2glu.2_only"),
  galactose = make_cond("scgal.1/fastlmm/results", "scgal.1"))
conds <- names(CONDS)

# ---- strong set: lead variant per gene (does not depend on the seed) ----
signif <- lapply(CONDS, function(cfg) bind_rows(lapply(VARIANTS, function(vt)
  fread(cfg$signif[[vt]], sep = "\t", quote = "", data.table = FALSE) %>%
    transmute(Gene = Source_Phenotype, SNP, Variant_Type = vt, P = PValue))) %>%
  distinct(Gene, SNP, Variant_Type, .keep_all = TRUE))
strong <- bind_rows(signif) %>% group_by(Gene) %>% slice_min(P, n = 1, with_ties = FALSE) %>% ungroup() %>%
  select(Gene, SNP, Variant_Type)

# ---- random set: the reproducible draw of the main analysis, with this seed ----
N_RANDOM_GENES <- 1000; N_VAR_PER_GENE <- 5
random_pairs <- bind_rows(lapply(seq_along(VARIANTS), function(v) {
  vt <- VARIANTS[v]; d <- CONDS$YPD.1$dirs[[vt]]
  genes <- sort(list.dirs(d, full.names = FALSE, recursive = FALSE))
  set.seed(SEED + v)
  genes <- sample(genes, min(N_RANDOM_GENES, length(genes)))
  bind_rows(mclapply(seq_along(genes), function(i) {
    g <- genes[i]
    f <- file.path(d, g, paste0(g, ".first_assoc.txt"))
    if (!file.exists(f)) return(NULL)
    snps <- tryCatch(sort(fread(f, select = "SNP", showProgress = FALSE)$SNP), error = function(e) NULL)
    if (!length(snps)) return(NULL)
    set.seed(SEED + 1e5 * v + i)
    tibble(Gene = g, SNP = sample(snps, min(N_VAR_PER_GENE, length(snps))), Variant_Type = vt)
  }, mc.cores = N_CORES))
}))

# ---- effect sizes: reuse the main analysis' extraction, read only the pairs it does not have ----
key <- function(d) paste(d$Gene, d$SNP, d$Variant_Type, sep = "|")
need <- bind_rows(strong, random_pairs) %>% distinct()
extracted <- lapply(setNames(conds, conds), function(cn) {
  base <- readRDS(file.path(MASH_DIR, paste0("extract_", cn, ".rds")))
  todo <- need[!key(need) %in% key(base), ]
  jobs <- todo %>% group_by(Variant_Type, Gene) %>% summarise(snps = list(SNP), .groups = "drop")
  add <- rbindlist(mclapply(seq_len(nrow(jobs)), function(i) {
    vt <- jobs$Variant_Type[i]; g <- jobs$Gene[i]
    f <- file.path(CONDS[[cn]]$dirs[[vt]], g, paste0(g, ".first_assoc.txt"))
    if (!file.exists(f)) return(NULL)
    d <- tryCatch(fread(f, select = c("SNP", "SnpWeight", "SnpWeightSE", "PValue"), showProgress = FALSE), error = function(e) NULL)
    if (is.null(d)) return(NULL)
    d <- d[SNP %in% jobs$snps[[i]]]; if (!nrow(d)) return(NULL)
    d[, `:=`(Gene = g, Variant_Type = vt)]
  }, mc.cores = N_CORES))
  bind_rows(base, as.data.frame(add)) %>% distinct(Gene, SNP, Variant_Type, .keep_all = TRUE)
})

# ---- fit: identical to fit_run() of the main analysis ----
make_mats <- function(pairs) {
  k <- key(pairs)
  B <- S <- matrix(NA_real_, length(k), length(conds), dimnames = list(k, conds))
  for (cn in conds) { e <- extracted[[cn]]; m <- match(k, key(e)); B[, cn] <- e$SnpWeight[m]; S[, cn] <- e$SnpWeightSE[m] }
  list(B = B, S = S, bad = is.na(B) | is.na(S) | !is.finite(S) | S <= 0)
}
ms <- make_mats(strong); ms$B[ms$bad] <- 0; ms$S[ms$bad] <- 1e3
mr <- make_mats(random_pairs); keep <- rowSums(mr$bad) == 0
mr$B <- mr$B[keep, , drop = FALSE]; mr$S <- mr$S[keep, , drop = FALSE]
Vhat     <- estimate_null_correlation_simple(mash_set_data(mr$B, mr$S))
d_random <- mash_set_data(mr$B, mr$S, V = Vhat)
d_strong <- mash_set_data(ms$B, ms$S, V = Vhat)
U_ed <- cov_ed(d_strong, cov_pca(d_strong, npc = min(5, length(conds) - 1)))
m_random <- mash(d_random, Ulist = c(U_ed, cov_canonical(d_random)), outputlevel = 1)
m_strong <- mash(d_strong, g = get_fitted_g(m_random), fixg = TRUE)

saveRDS(list(seed = SEED, conds = conds, n_random_drawn = nrow(random_pairs), n_random = nrow(mr$B),
             pi = get_estimated_pi(m_random), Ulist = m_random$fitted_g$Ulist,
             pm = get_pm(m_strong), lfsr = get_lfsr(m_strong)),
        file.path(OUT_DIR, sprintf("seed_%d.rds", SEED)))
cat(sprintf("seed %d: %d random pairs used, null weight %.2f%%\n", SEED, nrow(mr$B), 100 * get_estimated_pi(m_random)[["null"]]))
