#!/bin/sh
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




# Build a per-matrix --keep list of the panel strains, taking the FID/IID pair
# straight from each matrix's own .fam so the format matches regardless of
# whether a strain is stored as a combined IID (SACE_YCR) or split FID+IID
# (SACE  YCR). Source of strain names = sample_used.txt (the panel definition).
# NOTE: LD/clumping only needs the panel strain LIST -- not the phenotype file.



set -euo pipefail

SAMPLES=/projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/sample_used.txt  # one strain name per line
OUTDIR=keep_lists
mkdir -p "$OUTDIR"

[ -f "$SAMPLES" ] || { echo "ERROR: $SAMPLES not found"; exit 1; }
n_samples=$(grep -c . "$SAMPLES")
echo "Panel definition: $n_samples strain names in $SAMPLES"

# One .fam per variant type
declare -A FAM=(
  [SNPs]=/projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/pan_transcriptome_paper/SNPs/SNPs_no_change_name_and_no_strain_filtering_only_MAF15_filtered/full969Matrix.SNPs.Biallelic.Var.DP10.99percNonMissing.MAF15.NoSubTelo20kb.plink.fam
  [CNVs]=/projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/pan_transcriptome_paper/CNVs/CNVs_no_change_name_and_no_strain_filtering_only_MAF15_filtered/CNV_969strains_gain_maf15_N400_20221121.plink.fam
  [Indels]=/projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/vcf_files_sent_by_Anne/SNPs_indels/modified_snp_indel_plink_for_96strains/Indels_only_96/indels.96samples.fam
  [SVs]=/projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/vcf_files_sent_by_Anne/SVs_from_SV.1087/SV.1087samples_updated_fam_files/SV.1087samples.SVs.fam
)

echo -e "\nType\tfam_total\tkeep_matched\tkeep_file"
for t in SNPs CNVs Indels SVs; do
  f="${FAM[$t]}"
  if [ ! -f "$f" ]; then echo -e "${t}\tFAM NOT FOUND: $f"; continue; fi
  keep="$OUTDIR/keep_${t}.txt"
  # Pull FID/IID straight from the .fam for rows whose FID, IID, or FID_IID
  # matches a panel strain name. Preserves that matrix's own FID/IID format.
  awk 'NR==FNR{ n[$1]=1; next }
       ($1 in n)||($2 in n)||(($1"_"$2) in n){ print $1"\t"$2 }' \
       "$SAMPLES" "$f" | sort -u > "$keep"
  fam_n=$(grep -c . "$f")
  keep_n=$(grep -c . "$keep")
  echo -e "${t}\t${fam_n}\t${keep_n}\t${keep}"
done

echo -e "\n'keep_matched' is the panel strains actually present in each matrix."
echo "Feed  keep_lists/keep_<TYPE>.txt  to BOTH  --keep  in run_plink_ld.sh"
echo "and to the matching  plink --clump  command for that type."
echo "If a type matches < $n_samples, those strains simply aren't genotyped there."
