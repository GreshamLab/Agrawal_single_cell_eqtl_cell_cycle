#!/usr/bin/env python3
# Usage: python3 generate_condition_rmds.py conditions.tsv
# conditions.tsv columns: Figure.4 folder name, tag in all_signif_* file names, FastLMM results dir,
# gene_location.txt, Seurat .rds (paths relative to common_analysis_for_all)
"""Create Figure.4/Fig.4-<condition>/{codes,results} with the three YPD.1 Rmds re-pathed per condition.

2026-10-16: the templates now use per-gene LD clumping in R (shared_code/ld_clump_and_local.R;
r2 > 0.6 within 50 kb) instead of PLINK, and classify eQTLs as local (25 kb; CNVs same ORF) or
distant. Panel A has a cis/trans map and a local/distant map; Step3 verifies local eQTLs."""
import os, re, sys

BASE = "/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all"
F3 = f"{BASE}/Agrawal_single_cell_eqtl_cell_cycle/Figure.3"
F4 = f"{BASE}/Agrawal_single_cell_eqtl_cell_cycle/Figure.4"
F3_RESULTS = f"{F3}/Fig.3-Panel-A-C-D-cis-and-trans-eqtl-YPD-variant-type/results"

TEMPLATES = {
    f"{F3}/Fig.3-Panel-A-C-D-cis-and-trans-eqtl-YPD-variant-type/codes/Fig3_panel_B_C_D_variant_eqtl_classification.Rmd":
        "Step1_Fig4_panel_B_C_D_variant_eqtl_classification_{d}.Rmd",
    f"{F3}/Fig.3-Panel-A-C-D-cis-and-trans-eqtl-YPD-variant-type/codes/Fig.3-panel-A-cis_trans_eqtl_mapping_YPD.Rmd":
        "Step2_Fig4-panel-A-cis_trans_eqtl_mapping_{d}.Rmd",
    f"{F3}/Fig.3-Panel-B-raw_data_proof_YPD/Fig3-Panel-B-YPD.1_raw_eqtl_verification.Rmd":
        "Step3_Fig4-Panel-B-{d}_raw_eqtl_verification.Rmd",
}

def convert(text, d, tag, res, geneloc, seurat):
    out_results = f"{F4}/Fig.4-{d}/results"
    subs = [
        # output locations (Figure.3 YPD.1 results -> Figure.4 condition results)
        (f"{F3_RESULTS}/YPD.1", f"{out_results}/{d}"),
        (F3_RESULTS, out_results),
        # inputs: FastLMM results, gene locations, Seurat object
        (f"{BASE}/YPD.1/fastlmm/results", f"{BASE}/{res}"),
        (f"{BASE}/YPD.1/fastlmm/data/pheno_data/gene_location.txt", f"{BASE}/{geneloc}"),
        (f"{BASE}/YPD.1/gene_expression_analysis/output_from_pseudobulk_fastlmm/Na_va_data_subsetted.rds",
         f"{BASE}/{seurat}"),
        # association file names; non-YPD.1 InDel runs carry the _96_MAF15 tag
        ("all_signif_Indels_YPD.1_with_double_id_final_filtered",
         f"all_signif_Indels_{tag}_96_MAF15_with_double_id_final_filtered"),
        ("_YPD.1_", f"_{tag}_"),
        ('CONDITION <- "YPD.1"', f'CONDITION <- "{d}"'),
        ("all_clumped_eqtl_annotated_YPD.1.txt", f"all_clumped_eqtl_annotated_{d}.txt"),
        # titles
        ('title: "Fig.3-panel-A-cis_trans_eqtl_mapping_YPD"', f'title: "Fig.4-panel-A-cis_trans_eqtl_mapping_{d}"'),
        ('title: "Fig3-step1-_eQTL_Plotting_Analysis"', f'title: "Fig4-step1-eQTL_Plotting_Analysis_{d}"'),
        ('title: "Panel B — Single-cell raw-data validation of local eQTLs (all variant types)"',
         f'title: "{d} — Single-cell raw-data validation of local eQTLs (all variant types)"'),
        ("Fig3_panelA_cis_trans_2x2", f"Fig4_{d}_panelA_cis_trans_2x2"),
        ("Fig3_panelA_local_distant_2x2", f"Fig4_{d}_panelA_local_distant_2x2"),
        # Step1 saves annotated tables to PREV_ANALYSIS/annotations: make sure it exists
        ("for (dir in c(PREV_ANALYSIS, OUT_PLINK,", 'for (dir in c(PREV_ANALYSIS, file.path(PREV_ANALYSIS, "annotations"), OUT_PLINK,'),
        # a variant type with no local (or no distant) eQTLs has no local-vs-distant test
        ('summarise(p = wilcox.test(abs_SnpWeight ~ Location)$p.value, .groups = "drop")',
         'summarise(p = if (n_distinct(Location) == 2) wilcox.test(abs_SnpWeight ~ Location)$p.value else NA_real_,\n            .groups = "drop")'),
        # a condition can have fewer than 3 tested local eQTLs: no correlation then
        ('ct_p <- cor.test(pv$FastLMM_NegLog10P, pv$Verification_NegLog10P, method = "spearman", exact = FALSE)\n'
         'cat("Spearman rho (-log10 P, FastLMM vs NB GLM) =", round(ct_p$estimate, 3),\n'
         '    ", P =", signif(ct_p$p.value, 3), "\\n")',
         'if (nrow(pv) >= 3) {\n'
         '  ct_p <- cor.test(pv$FastLMM_NegLog10P, pv$Verification_NegLog10P, method = "spearman", exact = FALSE)\n'
         '  cat("Spearman rho (-log10 P, FastLMM vs NB GLM) =", round(ct_p$estimate, 3),\n'
         '      ", P =", signif(ct_p$p.value, 3), "\\n")\n'
         '} else cat("Spearman rho not computed: only", nrow(pv), "tested local eQTL(s)\\n")'),
        # integrated (replicate-merged) objects keep one counts layer per replicate
        ("seurat_object <- readRDS(SEURAT_FILE)\n",
         "seurat_object <- readRDS(SEURAT_FILE)\n"
         "# integrated objects store counts per replicate (counts.1, counts.2): join them\n"
         "if (!\"counts\" %in% Layers(seurat_object[[\"RNA\"]]))\n"
         "  seurat_object[[\"RNA\"]] <- JoinLayers(seurat_object[[\"RNA\"]])\n"),
    ]
    for a, b in subs:
        text = text.replace(a, b)
    return text

def main(table):
    for line in open(table):
        if not line.strip():
            continue
        d, tag, res, geneloc, seurat = line.rstrip("\n").split("\t")
        root = f"{F4}/Fig.4-{d}"
        for sub in ("codes", "results"):
            os.makedirs(f"{root}/{sub}", exist_ok=True)
        for src, dst in TEMPLATES.items():
            text = convert(open(src).read(), d, tag, res, geneloc, seurat)
            left = [l for l in text.splitlines() if re.search(r"YPD\.1|Figure\.3", l)]
            path = f"{root}/codes/{dst.format(d=d)}"
            with open(path, "w") as fh:
                fh.write(text)
            print(f"wrote {path}" + (f"  [WARN {len(left)} YPD.1/Figure.3 refs left]" if left else ""))
            for l in left:
                print("    ", l.strip()[:160])

if __name__ == "__main__":
    main(sys.argv[1])
