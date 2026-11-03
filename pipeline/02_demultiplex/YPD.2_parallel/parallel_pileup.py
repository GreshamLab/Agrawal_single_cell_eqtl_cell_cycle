import subprocess
import os
import sys
from multiprocessing import Pool

index = int(sys.argv[1])
cell_barcodes_path = f"../split_bam_new_code/barcodes/cell_barcodes_{index}.txt"
popscle_sif = "/projects/rps/cgsb/gresham/Akriti/wildyeast1/moved_from_scratch/code/popscle_2021-05-05-da70fc7.sif"
vcf_file = "/projects/rps/cgsb/gresham/Akriti/wildyeast1/moved_from_scratch/less_stringent_demuxlet/filteredtill_mutations_homo_in_all.vcf"
bam_dir = "../split_bam_new_code/bams_YPD.2_new"
output_base_dir = "./pileup_YPD.2"
os.makedirs(output_base_dir, exist_ok=True)

with open(cell_barcodes_path, 'r') as f:
    barcodes = [line.strip() for line in f if line.strip()]

def run_pileup(barcode):
    bam_filename = f"{bam_dir}/{barcode}.bam"
    barcode_output_dir = f"{output_base_dir}/{barcode}"
    os.makedirs(barcode_output_dir, exist_ok=True)
    output_pileup = f"{barcode_output_dir}/samples_to_demuxlet.pileup"
    command = f"singularity exec {popscle_sif} popscle dsc-pileup --sam {bam_filename} --vcf {vcf_file} --group-list {cell_barcodes_path} --out {output_pileup}"
    subprocess.run(command, shell=True)

# Run N barcodes in parallel within each job
# Use same number as --cpus-per-task in your SLURM script
with Pool(processes=3) as pool:
    pool.map(run_pileup, barcodes)
