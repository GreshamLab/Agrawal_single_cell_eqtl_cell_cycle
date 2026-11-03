#!/bin/bash
# Gene-CNV genotypes without chromosome 17 (unplaced CNV markers), used for the gene-CNV FaST-LMM
# re-run of 2026-10-22. Command as recorded in CNVs.96strains.maf15.no_chr17.plink.log (PLINK v1.9).
CNV=/projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/structural_variant_vcf/SVcalls2/SVcalls/CNVs_96_strains/CNV_96strains_maf15
mkdir -p $CNV/no_chr17 && cd $CNV/no_chr17
plink --bfile $CNV/CNVs.96strains.maf15.plink --chr 1-16 --make-bed --out CNVs.96strains.maf15.no_chr17.plink
