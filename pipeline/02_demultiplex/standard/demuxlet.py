import subprocess
import os
import sys

# Get the index from the command line argument
index = int(sys.argv[1])

# Directory where barcodes are stored
cell_barcodes_path = f"../split_bam/barcodes/cell_barcodes_{index}.txt"

# Full barcode list path
full_barcode_list = "../gene_expression_analysis/cell_barcodes_1SD.txt"

# Constants for the popscle demuxlet command
popscle_sif = "/scratch/cgsb/gresham/Akriti/wildyeast1/code/popscle_2021-05-05-da70fc7.sif"
vcf_file = "/scratch/aa9679/SingleCell_cerevisiae/data.2/less_stringent_demuxlet/filteredtill_mutations_homo_in_all.vcf"

# Directory where pileup files are stored
pileup_base_dir = "../pileup/pileup_YPD1_out_20k_1SD"

# New directory for demuxlet output
output_base_dir = "./demuxlet_YPD1_out_1SD_default"
os.makedirs(output_base_dir, exist_ok=True)

# Ensure the VCF file exists
if not os.path.exists(vcf_file):
    sys.exit(f"VCF file not found: {vcf_file}")

# Ensure the full barcode list file exists
if not os.path.exists(full_barcode_list):
    sys.exit(f"Full barcode list file not found: {full_barcode_list}")

# Open and process each barcode in the file
try:
    with open(cell_barcodes_path, 'r') as barcode_file:
        for line in barcode_file:
            barcode = line.strip()

            # Construct the path to the required pileup file
            plp_file = f"{pileup_base_dir}/{barcode}/samples_to_demuxlet.pileup"

            # Print the path for debugging
            print(f"Checking for pileup file: {plp_file}")

            # Verify that the required pileup file exists
            if not os.path.exists(plp_file + '.plp.gz'):
                print(f"Missing pileup file for barcode: {barcode}")
                continue


            # Create a directory to store demuxlet output for this barcode in the new output directory
            barcode_output_dir = f"{output_base_dir}/{barcode}"
            os.makedirs(barcode_output_dir, exist_ok=True)

            # Construct the output demuxlet TSV file path in the new output directory
            output_demuxlet = f"{barcode_output_dir}/demuxlet_output.tsv"

            # Singularity popscle demuxlet command with the plp file only
            demuxlet_command = (
                f"singularity exec {popscle_sif} popscle demuxlet "
                f"--plp {plp_file} --field GT --vcf {vcf_file} "
                f"--group-list {full_barcode_list} --out {output_demuxlet} "
            )

            # Print the command for debugging
            print(f"Running demuxlet command: {demuxlet_command}")

            # Run the demuxlet command
            result = subprocess.run(demuxlet_command, shell=True, capture_output=True, text=True)

            # Check if the command was successful
            if result.returncode != 0:
                print(f"Error running demuxlet for barcode {barcode}:\n{result.stderr}")
            else:
                print(f"Demuxlet completed successfully for barcode {barcode}")

except FileNotFoundError as e:
    print(f"Error: {e}")


