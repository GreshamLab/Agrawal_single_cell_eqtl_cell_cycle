#!/bin/bash
# Re-run every Figure 3 / Figure 4 eQTL analysis downstream of FaST-LMM.
# FaST-LMM itself and the cached PLINK clumps / mashr fits are reused (inputs unchanged).
# Run from this directory:  bash submit_all.sh
cd "$(dirname "$0")"
ids=()
for ds in YPD.1 YPD.2 YPD_integrated Flu.1 Flu.2 Flu_integrated sc2glu.1_only sc2glu.2_only sc2glu_integrated scgal.1; do
  id=$(sbatch --parsable -J "eqtl_${ds}" run_dataset.sbatch "$ds")
  echo "$ds -> $id"; ids+=("$id")
done
dep=$(IFS=:; echo "${ids[*]}")
id=$(sbatch --parsable -J eqtl_cross_dataset --dependency=afterany:$dep run_cross_dataset.sbatch)
echo "cross-dataset -> $id (starts after $dep)"
