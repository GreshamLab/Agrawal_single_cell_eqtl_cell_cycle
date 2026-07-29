#!/bin/bash
#
#SBATCH --verbose
#SBATCH --job-name=LD
#SBATCH --output=LD.o%j
#SBATCH --error=LD.e%j
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=10gb
#SBATCH --mail-type=END
#SBATCH --account=torch_pr_130_general
#SBATCH --mail-user=aa9679@nyu.edu

module load bioinformatics/20260224

# Emit one pairwise LD table (.ld) per variant type, on the panel strains only.
# Run build_keep_and_check.sh FIRST to create keep_lists/keep_<TYPE>.txt.
set -uo pipefail   # no -e: one failed type should not abort the rest

WINDOW_KB=20       # max distance to measure (go past where r2 flattens)

KEEPDIR=keep_lists
if [ ! -d "$KEEPDIR" ]; then
  echo "ERROR: $KEEPDIR/ not found - run build_keep_and_check.sh first" >&2
  exit 1
fi

mkdir -p ld_tables

while read -r NAME BFILE; do
  [ -z "$NAME" ] && continue
  KEEP="$KEEPDIR/keep_${NAME}.txt"
  if [ ! -f "$KEEP" ]; then echo "SKIP $NAME: no $KEEP" >&2; continue; fi
  echo "=== $NAME  (keep $(grep -c . "$KEEP") strains) ==="
  plink --bfile "$BFILE" \
        --keep "$KEEP" \
        --r2 \
        --ld-window-kb  "$WINDOW_KB" \
        --ld-window     999999 \
        --ld-window-r2  0 \
        --out ld_tables/ld_${NAME} </dev/null
  if [ -f "ld_tables/ld_${NAME}.ld" ]; then
    echo "  OK -> ld_tables/ld_${NAME}.ld ($(wc -l < ld_tables/ld_${NAME}.ld) lines)"
  else
    echo "  FAILED for $NAME - see ld_tables/ld_${NAME}.log" >&2
  fi
done <<'EOF'
SNPs    /projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/pan_transcriptome_paper/SNPs/SNPs_no_change_name_and_no_strain_filtering_only_MAF15_filtered/full969Matrix.SNPs.Biallelic.Var.DP10.99percNonMissing.MAF15.NoSubTelo20kb.plink
CNVs    /projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/pan_transcriptome_paper/CNVs/CNVs_no_change_name_and_no_strain_filtering_only_MAF15_filtered/CNV_969strains_gain_maf15_N400_20221121.plink
Indels  /projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/vcf_files_sent_by_Anne/SNPs_indels/modified_snp_indel_plink_for_96strains/Indels_only_96/indels.96samples
SVs     /projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/vcf_files_sent_by_Anne/SVs_from_SV.1087/SV.1087samples_updated_fam_files/SV.1087samples.SVs
EOF

echo "Done. .ld tables in ./ld_tables/"
echo "REMINDER: add  --keep keep_lists/keep_<TYPE>.txt  to each plink --clump too."