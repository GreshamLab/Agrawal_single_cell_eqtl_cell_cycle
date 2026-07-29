#!/usr/bin/env python3
import os
import glob
import math

def create_phenotype_list(pheno_dir, output_file):
    """Create phenotype file list"""
    
    pheno_dir = os.path.abspath(pheno_dir)
    print(f"Searching directory: {pheno_dir}")
    
    # Get all .txt files and sort them
    pheno_files = sorted(glob.glob(os.path.join(pheno_dir, "*.norm.txt")))
    
    print(f"Found {len(pheno_files)} phenotype files")
    
    if len(pheno_files) == 0:
        print("ERROR: No phenotype files found!")
        return None, 0
    
    # Write absolute paths to file
    with open(output_file, 'w') as f:
        for pheno_file in pheno_files:
            f.write(os.path.abspath(pheno_file) + '\n')
    
    print(f"Written to: {output_file}")
    print(f"First 3 files: {', '.join([os.path.basename(f) for f in pheno_files[:3]])}")
    print(f"Last 3 files: {', '.join([os.path.basename(f) for f in pheno_files[-3:]])}")
    
    return output_file, len(pheno_files)

def split_into_batches(list_file, lines_per_batch, output_dir):
    """Split phenotype list into batches"""
    
    print(f"\n--- Splitting into batches ---")
    
    # Read all lines
    with open(list_file, 'r') as f:
        lines = f.readlines()
    
    total_lines = len(lines)
    total_batches = math.ceil(total_lines / lines_per_batch)
    
    print(f"Total files: {total_lines}")
    print(f"Files per batch: {lines_per_batch}")
    print(f"Creating: {total_batches} batch files\n")
    
    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)
    
    # Split into batch files
    for i in range(total_batches):
        start_idx = i * lines_per_batch
        end_idx = min((i + 1) * lines_per_batch, total_lines)
        batch_lines = lines[start_idx:end_idx]
        
        batch_file = os.path.join(output_dir, f"phenotype_batch_{i+1:02d}.txt")
        
        with open(batch_file, 'w') as f:
            f.writelines(batch_lines)
        
        print(f"Batch {i+1:02d}: {len(batch_lines)} files → {os.path.basename(batch_file)}")
    
    return total_batches

if __name__ == "__main__":
    # Configuration
    pheno_dir = os.path.abspath("../data/normalized_pheno_file")
    list_output_dir = os.path.abspath("../data/phenotype_lists")
    batch_output_dir = os.path.abspath("../data/pheno_batches")
    lines_per_batch = 10  # Files per batch
    
    # Create output directories
    os.makedirs(list_output_dir, exist_ok=True)
    

    # Step 1: Create phenotype list
    list_file = os.path.join(list_output_dir, "phenotype_files.txt")
    list_file, total_files = create_phenotype_list(pheno_dir, list_file)
    
    if list_file is None:
        exit(1)
    
    # Step 2: Split into batches
    total_batches = split_into_batches(list_file, lines_per_batch, batch_output_dir)

