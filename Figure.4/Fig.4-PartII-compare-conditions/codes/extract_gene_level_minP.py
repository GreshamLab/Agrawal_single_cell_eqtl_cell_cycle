#!/usr/bin/env python3
"""
Gene-level eQTL statistics for every TESTED gene (not only significant ones), used as the
GO universe and for the ranked GO test in Fig4.PartII-compare-conditions-best-replicates.Rmd.

For each condition x variant type x gene, reads the FaST-LMM output directory used by the
Step1 classification Rmds:
  <gene>.first_assoc.txt  sorted by PValue, so row 1 is the gene's top variant (min p)
  <gene>.threshold.txt    5th percentile of 100 permuted min p-values (thr_perm)
min_p < thr_perm  <=>  gene-level permutation p < 0.05 for that variant type.

Output: ../results/best_replicates/gene_level_minP.tsv
Usage (from codes/): python3 extract_gene_level_minP.py
"""

import os
import csv

ROOT = "/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all"
SUF_NEW = "with_double_id_is_not_needed_final_filtered_for_ncount_less_than_5_run_from_new_plink_files"
SNP_96 = "results_SNPs_96_MAF15_from_969_strains_with_double_id_final_filtered_for_ncount_less_than_5"
IND_96 = "results_Indels_96_MAF15_with_double_id_final_filtered_for_ncount_less_than_5"

# same runs as the assoc files in each condition's Step1 classification Rmd
DIRS = {
    "YPD.1": ("YPD.1/fastlmm/results", {
        "SNPs": SNP_96,
        "InDels": "results_Indels_with_double_id_final_filtered_for_ncount_less_than_5",
        "CNVs": f"results_CNVs_{SUF_NEW}_no_chr17", "SVs": f"results_SVs_{SUF_NEW}"}),
    "sc2glu.2": ("sc2glu.2/fastlmm_sc2.glu.2_only/results", {
        "SNPs": SNP_96, "InDels": IND_96, "CNVs": f"results_CNVs_{SUF_NEW}_no_chr17", "SVs": f"results_SVs_{SUF_NEW}"}),
    "Flu.1": ("Flu.1/fastlmm/results", {
        "SNPs": SNP_96, "InDels": IND_96, "CNVs": f"results_CNVs_{SUF_NEW}_no_chr17", "SVs": f"results_SVs_{SUF_NEW}"}),
    "Flu.2": ("Flu.2/fastlmm/results", {
        "SNPs": SNP_96, "InDels": IND_96, "CNVs": f"results_CNVs_{SUF_NEW}_no_chr17", "SVs": f"results_SVs_{SUF_NEW}"}),
}

OUT = "../results/best_replicates/gene_level_minP.tsv"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["condition", "variant_type", "gene", "top_variant", "min_p", "effect", "thr_perm"])
    for cond, (base, vts) in DIRS.items():
        for vt, sub in vts.items():
            d = os.path.join(ROOT, base, sub)
            n = miss = 0
            for gene in sorted(os.listdir(d)):
                gdir = os.path.join(d, gene)
                if not os.path.isdir(gdir):
                    continue
                assoc = os.path.join(gdir, f"{gene}.first_assoc.txt")
                thr_f = os.path.join(gdir, f"{gene}.threshold.txt")
                if not (os.path.exists(assoc) and os.path.exists(thr_f)):
                    miss += 1
                    continue
                with open(assoc) as a:
                    header = a.readline().rstrip("\n").split("\t")
                    row = a.readline().rstrip("\n").split("\t")
                with open(thr_f) as t:
                    t.readline()
                    thr = t.readline().strip()
                if len(row) < len(header) or not thr:
                    miss += 1
                    continue
                r = dict(zip(header, row))
                w.writerow([cond, vt, gene, r["SNP"], r["PValue"], r["SnpWeight"], thr])
                n += 1
            print(f"{cond:9s} {vt:7s} genes={n:5d} missing/empty={miss}")
