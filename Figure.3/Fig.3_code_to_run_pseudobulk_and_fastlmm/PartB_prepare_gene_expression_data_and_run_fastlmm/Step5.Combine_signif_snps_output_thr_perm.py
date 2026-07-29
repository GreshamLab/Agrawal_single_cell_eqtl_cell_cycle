#!/usr/bin/env python3
"""
Add permutation threshold (thr_perm) values to significant SNPs data

Purpose: For each phenotype, read the threshold.txt file and add the 
thr_perm column to all significant SNPs for that phenotype

The threshold.txt file contains the p-value threshold from permutation testing
This threshold is used to determine if an eQTL is significant

Usage: python add_threshold_to_signif_snps.py
"""

import pandas as pd
import glob
import os

# Define variant types and their directories
variant_types = {
    'SNP_YPD.1_with_double_id': '../results/results_SNPs_968_with_double_id',
    'CNV_YPD.1_with_double_id': '../results/results_CNVs_with_double_id',
#    'SVs_84': '../results/results_SVs_84',
    'Indels_YPD.1_with_double_id': '../results/results_Indels_with_double_id',
}

# Process each variant type
for variant_name, results_dir in variant_types.items():
    print(f"\n{'=' * 80}")
    print(f"Processing {variant_name}")
    print(f"{'=' * 80}")
    
    # Find all signif_snps.txt files
    files = glob.glob(f"{results_dir}/*/*.signif_snps.txt")
    print(f"Found {len(files)} signif_snps.txt files\n")
    
    # Read and combine all files WITH THRESHOLD VALUES
    all_data = []
    empty_files = []
    phenotypes_with_missing_threshold = []
    
    for file in files:
        # Extract phenotype name from directory
        phenotype = os.path.basename(os.path.dirname(file))
        phenotype_dir = os.path.dirname(file)
        
        # Read the signif_snps file
        try:
            df = pd.read_csv(file, sep='\t')
            
            if not df.empty:
                # Read corresponding threshold file
                threshold_file = os.path.join(phenotype_dir, f"{phenotype}.threshold.txt")
                
                if os.path.exists(threshold_file):
                    try:
                        # Read threshold file
                        # Format is typically: x (header)\n<threshold_value>
                        threshold_df = pd.read_csv(threshold_file, sep='\t', header=0)
                        
                        # Extract threshold value (should be in first row, first column)
                        if len(threshold_df) > 0:
                            threshold_value = threshold_df.iloc[0, 0]
                            
                            # Validate threshold is a number
                            try:
                                threshold_value = float(threshold_value)
                            except ValueError:
                                print(f"  ⚠ Invalid threshold value for {phenotype}: {threshold_value}")
                                phenotypes_with_missing_threshold.append(phenotype)
                                continue
                            
                            # Add phenotype identifier and threshold to dataframe
                            df['Source_Phenotype'] = phenotype
                            df['thr_perm'] = threshold_value
                            
                            all_data.append(df)
                            print(f"  ✓ {phenotype:40s} | {len(df):5d} hits | thr_perm: {threshold_value:.2e}")
                        else:
                            print(f"  ✗ Empty threshold file for {phenotype}")
                            phenotypes_with_missing_threshold.append(phenotype)
                    except Exception as e:
                        print(f"  ✗ Error reading threshold for {phenotype}: {e}")
                        phenotypes_with_missing_threshold.append(phenotype)
                else:
                    print(f"  ✗ Threshold file not found: {threshold_file}")
                    phenotypes_with_missing_threshold.append(phenotype)
            else:
                empty_files.append(phenotype)
                print(f"  ⊘ Empty signif_snps file for {phenotype}")
                
        except Exception as e:
            print(f"  ✗ Could not read {file}: {e}")
    
    # Report issues
    print(f"\n{'─' * 80}")
    print(f"Summary of file processing:")
    print(f"{'─' * 80}")
    print(f"  ✓ Successfully processed: {len(all_data)} phenotypes with data")
    print(f"  ⊘ Empty signif_snps files: {len(empty_files)}")
    print(f"  ✗ Missing/error in threshold: {len(phenotypes_with_missing_threshold)}")
    
    if phenotypes_with_missing_threshold:
        print(f"\n  Phenotypes with threshold issues:")
        for pheno in phenotypes_with_missing_threshold[:10]:
            print(f"    - {pheno}")
        if len(phenotypes_with_missing_threshold) > 10:
            print(f"    ... and {len(phenotypes_with_missing_threshold) - 10} more")
    
    # Combine all data
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        
        # Reorder columns to put Source_Phenotype and thr_perm first
        cols = ['Source_Phenotype', 'thr_perm'] + [col for col in combined_df.columns 
                                                     if col not in ['Source_Phenotype', 'thr_perm']]
        combined_df = combined_df[cols]
        
        # Save combined file
        output_file = f"../results/all_signif_{variant_name}_with_threshold.txt"
        combined_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"\n{'=' * 80}")
        print(f"✓ Output saved to: {output_file}")
        print(f"{'=' * 80}")
        
        print(f"\nCombined {len(all_data)} files with {len(combined_df)} total significant hits")
        
        # Summary statistics
        print(f"\nSummary for {variant_name}:")
        print(f"  Total significant hits: {len(combined_df)}")
        print(f"  Phenotypes with significant hits: {combined_df['Source_Phenotype'].nunique()}")
        print(f"  Threshold p-value range: {combined_df['thr_perm'].min():.2e} to {combined_df['thr_perm'].max():.2e}")
        print(f"\n  Top 10 phenotypes by number of significant hits:")
        top_10 = combined_df['Source_Phenotype'].value_counts().head(10)
        for pheno, count in top_10.items():
            thr_val = combined_df[combined_df['Source_Phenotype'] == pheno]['thr_perm'].iloc[0]
            print(f"    {pheno:40s} | {count:5d} hits | thr_perm: {thr_val:.2e}")
        
    else:
        print(f"  No significant hits found for {variant_name}")

print(f"\n{'=' * 80}")
print("ALL VARIANT TYPES PROCESSED")
print(f"{'=' * 80}\n")
