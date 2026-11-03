# =============================================================================
# Pseudobulk count matrices WITHOUT the >5 cells per strain filter
#
# The >5 cells per strain filter is only for the eQTL step. The pseudobulk
# scripts of most datasets apply it before writing the pseudobulk matrix
# (pseudobulk_<ds>_without_filtering_transcripts.txt: "without filtering"
# refers to transcripts, not strains). This script re-runs each dataset's own
# pseudobulk Rmd (same cell QC and strain exclusions) with the cell threshold
# as a parameter, and stops once the raw pseudobulk matrix is written.
#
# All other writes in those Rmds (saveRDS, genotype lists, fastlmm files) are
# intercepted and skipped, so nothing used by the eQTL pipeline is touched.
#
# Usage: Rscript 0_pseudobulk_no_cell_filter.R <dataset> <min_cells> <out_file>
#   min_cells = 5 reproduces the existing (filtered) matrix, for validation
#   min_cells = 0 keeps every demultiplexed strain
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
DATASET   <- args[1]
MIN_CELLS <- as.numeric(args[2])
OUT_FILE  <- args[3]

BASE <- "/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all"

# dataset -> source Rmd and the file name its raw pseudobulk matrix is written to
SOURCES <- list(
  YPD.1   = c("YPD.1/gene_expression_analysis/Pseusdobulk_fastlmm.Rmd",       "pseudobulk_YPD.1_without_filtering_transcripts.txt"),
  YPD.2   = c("YPD.2/gene_expression_analysis/Pseusdobulk_fastlmm.Rmd",       "pseudobulk_YPD.2_without_filtering_transcripts.txt"),
  Flu.1   = c("Flu.1/gene_expression_analysis/Pseusdobulk_fastlmm_Flu1.Rmd",  "pseudobulk_flu.1_without_filtering_transcripts.txt"),
  Flu.2   = c("Flu.2/gene_expression_analysis/Pseusdobulk_fastlmm_flu2.Rmd",  "pseudobulk_flu.2_without_filtering_transcripts.txt"),
  scgal.1 = c("scgal.1/gene_expression_analysis/Pseusdobulk_fastlmm.Rmd",     "pseudobulk_scgal.1_without_filtering_transcripts.txt"),
  sc2glu_integrated = c("sc2glu.2/sc2glu.1_and_sc2glu.2_seurat_merged/merged_outputs_from_sc2glu.1_and_sc2glu.2.Rmd",
                        "pseudobulked_sc2glu_before_any_filtering.csv")
)
src    <- file.path(BASE, SOURCES[[DATASET]][1])
target <- SOURCES[[DATASET]][2]

# Rmd -> R code; make the cell threshold a parameter
r_file <- tempfile(fileext = ".R")
knitr::purl(src, output = r_file, quiet = TRUE, documentation = 0)
code <- readLines(r_file)
n_filter <- sum(grepl("n_cells > 5", code, fixed = TRUE))
cat(DATASET, ": replacing", n_filter, "'n_cells > 5' filter(s) with n_cells >", MIN_CELLS, "\n")
stopifnot(n_filter >= 1)
code <- gsub("n_cells > 5", "n_cells > MIN_CELLS", code, fixed = TRUE)
# the Rmds were run with SeuratObject < 5, where GetAssayData() took `slot`
code <- gsub("slot = \"counts\"", "layer = \"counts\"", code, fixed = TRUE)
exprs <- parse(text = code)

env <- new.env()
env$MIN_CELLS <- MIN_CELLS
env$done <- FALSE

# intercept writes: only the target pseudobulk matrix is written (to OUT_FILE)
intercept <- function(writer) {
  force(writer)
  function(x, file = "", ...) {
    if (is.character(file) && basename(file) == target) {
      writer(x, file = OUT_FILE, ...)
      cat("Wrote", ncol(x), "strains x", nrow(x), "genes to", OUT_FILE, "\n")
      env$done <- TRUE
    } else {
      cat("Skipped write:", if (is.character(file)) file else "<connection>", "\n")
    }
    invisible(NULL)
  }
}
env$write.table <- intercept(utils::write.table)
env$write.csv   <- intercept(utils::write.csv)
env$saveRDS     <- function(object, file = "", ...) { cat("Skipped saveRDS:", file, "\n"); invisible(NULL) }

setwd(dirname(src))   # the Rmds use paths relative to their own folder
for (e in exprs) {
  eval(e, envir = env)
  if (env$done) break
}
stopifnot(env$done)
