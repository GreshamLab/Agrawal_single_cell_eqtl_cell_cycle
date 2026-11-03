# ============================================================
# Summary of the mashr fits with five random sets (mashr_seed_robustness.R): pattern weights,
# data-driven patterns, sharing and condition-specific calls per seed, their range, and the
# eQTLs whose call is the same in most fits. Strong-set eQTLs without a local / distant call are
# dropped, as in the main analysis.
# Output: results/mashr/seed_robustness/
#   seed_summary.csv          one row per seed (and median / min / max)
#   seed_patterns.csv         every pattern with a non-zero weight, per seed
#   seed_consensus_calls.csv  number of eQTLs per call made in >= 3 of 5 fits
# ============================================================
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
select <- dplyr::select; filter <- dplyr::filter

PROJ <- "/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/Agrawal_single_cell_eqtl_cell_cycle"
OUT_DIR <- file.path(PROJ, "Figure.4/Fig.4-PartII-compare-conditions/results/mashr/seed_robustness")
SEEDS <- c(42, 101, 202, 303, 404)
fits <- lapply(setNames(SEEDS, SEEDS), function(s) readRDS(file.path(OUT_DIR, sprintf("seed_%d.rds", s))))

# seed 42 must reproduce the main fit
main <- readRDS(file.path(PROJ, "Figure.4/Fig.4-PartII-compare-conditions/results/mashr/mash_A_main.rds"))
pi_main <- mashr::get_estimated_pi(main$m_random)
cat("seed 42 reproduces the main fit:",
    isTRUE(all.equal(unname(pi_main), unname(fits[["42"]]$pi), tolerance = 1e-6)), "\n")

# local / distant calls (as in the main analysis), to drop eQTLs without a call
ANNOT_DIRS <- c(YPD.1 = "Figure.3/Fig.3-Panel-A-C-D-cis-and-trans-eqtl-YPD-variant-type/results/YPD.1/annotations",
                YPD.2 = "Figure.4/Fig.4-YPD.2/results/YPD.2/annotations",
                Flu_16ug = "Figure.4/Fig.4-Flu.1/results/Flu.1/annotations",
                Flu_8ug = "Figure.4/Fig.4-Flu.2/results/Flu.2/annotations",
                sc2glu.2 = "Figure.4/Fig.4-sc2glu.2_only/results/sc2glu.2_only/annotations",
                galactose = "Figure.4/Fig.4-scgal.1/results/scgal.1/annotations")
ANNOT_FILES <- c(SNPs = "snp", InDels = "indel", SVs = "sv", CNVs = "cnv")
ld_calls <- bind_rows(lapply(names(ANNOT_DIRS), function(cn) bind_rows(imap(ANNOT_FILES, function(f, vt)
  fread(file.path(PROJ, ANNOT_DIRS[[cn]], paste0(f, "_data_annotated_cis_trans.txt")), sep = "\t", quote = "",
        select = c("Pheno", "SNP", "clump_members", "local_distant"), data.table = FALSE) %>%
    mutate(Variant_Type = vt))))) %>%
  mutate(member = strsplit(coalesce(clump_members, SNP), ";")) %>% unnest(member) %>%
  distinct(Gene = Pheno, SNP = member, Variant_Type, local_distant) %>%
  mutate(key = paste(Gene, SNP, Variant_Type, sep = "|"))

kind_of <- function(p, conds) case_when(p == "null" ~ "no_effect", grepl("^ED_", p) ~ "data_driven",
                                        p %in% c("equal_effects", "identity") | grepl("^simple_het", p) ~ "equal_effects",
                                        p %in% conds ~ "one_dataset", TRUE ~ "other")

per_seed <- lapply(fits, function(f) {
  conds <- f$conds
  ok <- rownames(f$pm) %in% ld_calls$key
  pm <- f$pm[ok, , drop = FALSE]; lf <- f$lfsr[ok, , drop = FALSE]
  sig <- rowSums(lf < 0.05) > 0
  rel <- t(apply(pm[sig, , drop = FALSE], 1, function(x) x / x[which.max(abs(x))]))
  loc <- ld_calls$local_distant[match(rownames(pm), ld_calls$key)]
  sg <- function(c) lf[, c] < 0.05
  sh <- function(a, b) sign(pm[, a]) == sign(pm[, b]) & abs(pm[, b] / pm[, a]) >= 0.5 & abs(pm[, b] / pm[, a]) <= 2
  nonY <- setdiff(conds, c("YPD.1", "YPD.2")); nonF <- setdiff(conds, c("Flu_16ug", "Flu_8ug"))
  calls <- list(
    YPD_specific = sg("YPD.1") & sg("YPD.2") & Reduce(`&`, lapply(nonY, function(o) !sh("YPD.1", o) & !sh("YPD.2", o))),
    fluconazole_specific = sg("Flu_16ug") & sg("Flu_8ug") & Reduce(`&`, lapply(nonF, function(o) !sh("Flu_16ug", o) & !sh("Flu_8ug", o))),
    galactose_specific = sg("galactose") & Reduce(`&`, lapply(setdiff(conds, "galactose"), function(o) !sh("galactose", o))),
    SC_glucose_specific = sg("sc2glu.2") & Reduce(`&`, lapply(setdiff(conds, "sc2glu.2"), function(o) !sh("sc2glu.2", o))),
    noise_YPD1_not_YPD2 = sg("YPD.1") & !sh("YPD.1", "YPD.2"),
    noise_YPD2_not_YPD1 = sg("YPD.2") & !sh("YPD.2", "YPD.1"))
  calls <- lapply(calls, function(x) x & !is.na(x))
  mag_all6 <- setNames(rep(FALSE, nrow(pm)), rownames(pm)); mag_all6[rownames(rel)] <- rowSums(rel >= 0.5) == length(conds)
  w <- tibble(pattern = names(f$pi), weight = 100 * f$pi) %>% mutate(kind = kind_of(pattern, conds))
  ed <- w %>% filter(kind == "data_driven", weight > 0) %>% arrange(desc(weight))
  top <- ed$pattern[1]; cm <- cov2cor(f$Ulist[[top]]); off <- cm[upper.tri(cm)]
  list(
    summary = tibble(seed = f$seed, random_pairs = f$n_random,
                     pct_no_effect = sum(w$weight[w$kind == "no_effect"]),
                     pct_data_driven = sum(w$weight[w$kind == "data_driven"]),
                     pct_equal_effects = sum(w$weight[w$kind == "equal_effects"]),
                     pct_one_dataset = sum(w$weight[w$kind == "one_dataset"]),
                     pct_shared_total = pct_data_driven + pct_equal_effects,
                     n_data_driven_patterns = nrow(ed), top_pattern = top, top_pattern_weight = ed$weight[1],
                     top_pattern_r_min = min(off), top_pattern_r_max = max(off),
                     strong_set = nrow(pm), confirmed = sum(sig),
                     local = sum(sig & loc == "Local", na.rm = TRUE), distant = sum(sig & loc == "Distant", na.rm = TRUE),
                     same_sign_all6 = sum(rowSums(rel > 0) == length(conds)),
                     pct_same_sign_all6 = 100 * mean(rowSums(rel > 0) == length(conds)),
                     sign_one_dataset = sum(rowSums(rel > 0) == 1),
                     magnitude_all6 = sum(rowSums(rel >= 0.5) == length(conds)),
                     pct_magnitude_all6 = 100 * mean(rowSums(rel >= 0.5) == length(conds)),
                     !!!lapply(calls, sum)),
    patterns = w %>% filter(weight > 0) %>% mutate(seed = f$seed),
    per_eqtl = tibble(eQTL = rownames(pm), confirmed = sig, magnitude_all6 = mag_all6, !!!calls))
})

summ <- bind_rows(lapply(per_seed, `[[`, "summary"))
num <- summ %>% select(where(is.numeric), -seed)
summ_out <- bind_rows(summ %>% mutate(seed = as.character(seed)),
                      tibble(seed = c("median", "min", "max")) %>%
                        bind_cols(bind_rows(num %>% summarise(across(everything(), median)),
                                            num %>% summarise(across(everything(), min)),
                                            num %>% summarise(across(everything(), max)))))
write.csv(summ_out, file.path(OUT_DIR, "seed_summary.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(per_seed, `[[`, "patterns")), file.path(OUT_DIR, "seed_patterns.csv"), row.names = FALSE)

# consensus: per eQTL, the number of fits (of 5) in which each call is made
pe <- bind_rows(lapply(per_seed, `[[`, "per_eqtl"), .id = "seed")
cons <- pe %>% select(-seed) %>% group_by(eQTL) %>% summarise(across(everything(), sum), .groups = "drop")
cons_counts <- cons %>% select(-eQTL) %>%
  summarise(across(everything(), list(in_all_5 = ~ sum(.x == 5), in_3plus = ~ sum(.x >= 3), in_any = ~ sum(.x >= 1)))) %>%
  pivot_longer(everything(), names_to = c("call", "how"), names_pattern = "(.*)_(in_.*)") %>%
  pivot_wider(names_from = how, values_from = value)
write.csv(cons_counts, file.path(OUT_DIR, "seed_consensus_calls.csv"), row.names = FALSE)
write.csv(cons, file.path(OUT_DIR, "seed_consensus_per_eQTL.csv"), row.names = FALSE)
print(summ_out, width = 250); print(cons_counts)
