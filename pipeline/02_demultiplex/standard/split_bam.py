ximport subprocess
import os
import sys

index = int(sys.argv[1])

# Path to the 'cell_barcodes' file
cell_barcodes_path = f"barcodes_Q4/cell_barcodes_{index}.txt"

#output directory for bam file
output_directory = "bams_Q4"
os.makedirs(output_directory, exist_ok=True)

# Iterate through each line in the file
with open(cell_barcodes_path, 'r') as file:
    for line in file:
        barcode = line.strip()
        
        # Create a temporary file for the current barcode
        temp_filename = f"temp_barcode_{index}_{barcode}.txt"
        with open(temp_filename, 'w') as temp_file:
            temp_file.write(barcode)
        
        # Construct the output BAM file name
        output_bam = os.path.join(output_directory, f"{barcode}.bam")
        
        # Run the subset-bam command
        command = f"/scratch/cgsb/gencore/mk5636/yeast_strain_analysis/subset-bam_linux --bam ../alignment/output_of_alignment_YPD1/run_count/outs/possorted_genome_bam.bam --cell-barcodes {temp_filename} --out-bam {output_bam} --cores 20"
        subprocess.run(command, shell=True)
        
        # Remove the temporary file
        os.remove(temp_filename)

