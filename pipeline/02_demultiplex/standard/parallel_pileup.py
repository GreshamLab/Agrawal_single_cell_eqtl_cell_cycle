import subprocess
import os
import sys

# Get the index from command line argument
index = int(sys.argv[1])

# Directory where barcodes are stored
# Path to the 'cell_barcodes' file, which will also be used as the group-list
cell_barcodes_path = f"../split_bam/barcodes_noSD/cell_barcodes_{index}.txt"

# Constants for the popscle command
popscle_sif = "/scratch/cgsb/gresham/Akriti/wildyeast1/code/popscle_2021-05-05-da70fc7.sif"
vcf_file = "/scratch/aa9679/SingleCell_cerevisiae/data.2/less_stringent_demuxlet/filteredtill_mutations_homo_in_all.vcf"

# Directory where BAM files are stored
bam_dir = "../split_bam/bams_1SD"

# Output directory for pileup files
output_base_dir = "./pileup_YPD1_out_noSD"
os.makedirs(output_base_dir, exist_ok=True)

# Open and process each barcode in the file
try:
    with open(cell_barcodes_path, 'r') as barcode_file:
        for line in barcode_file:
            barcode = line.strip()

            # Path to the associated BAM file for the barcode
            bam_filename = f"{bam_dir}/{barcode}.bam"

            # Create a subdirectory for each barcode inside the output directory
            barcode_output_dir = f"{output_base_dir}/{barcode}"
            os.makedirs(barcode_output_dir, exist_ok=True)

            # Construct the output pileup file path
            output_pileup = f"{barcode_output_dir}/samples_to_demuxlet.pileup"

            # Singularity popscle dsc-pileup command
            command = f"singularity exec {popscle_sif} popscle dsc-pileup --sam {bam_filename} --vcf {vcf_file} --group-list {cell_barcodes_path} --out {output_pileup}"

            # Run the command
            subprocess.run(command, shell=True)

except FileNotFoundError:
    print(f"Barcode file {cell_barcodes_path} not found.")
