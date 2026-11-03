import subprocess
import os
import sys
from multiprocessing import Pool

index = int(sys.argv[1])
cell_barcodes_path = f"../split_bam_new_code/barcodes/cell_barcodes_{index}.txt"
full_barcode_list = "../alignment/output_of_alignment_YPD2_with_clean_gtf/run_count/outs/filtered_feature_bc_matrix/barcodes.tsv"
popscle_sif = "/projects/rps/cgsb/gresham/Akriti/wildyeast1/moved_from_scratch/code/popscle_2021-05-05-da70fc7.sif"
vcf_file = "/projects/rps/cgsb/gresham/Akriti/wildyeast1/moved_from_scratch/less_stringent_demuxlet/filteredtill_mutations_homo_in_all.vcf"
pileup_base_dir = "../pileup/pileup_YPD.2"
output_base_dir = "./demuxlet_out_YPD.2"
os.makedirs(output_base_dir, exist_ok=True)

if not os.path.exists(vcf_file):
    sys.exit(f"VCF file not found: {vcf_file}")
if not os.path.exists(full_barcode_list):
    sys.exit(f"Full barcode list file not found: {full_barcode_list}")

with open(cell_barcodes_path, 'r') as f:
    barcodes = [line.strip() for line in f if line.strip()]

def run_demuxlet(barcode):
    plp_file = f"{pileup_base_dir}/{barcode}/samples_to_demuxlet.pileup"
    print(f"Checking for pileup file: {plp_file}")

    if not os.path.exists(plp_file + '.plp.gz'):
        print(f"Missing pileup file for barcode: {barcode}")
        return

    barcode_output_dir = f"{output_base_dir}/{barcode}"
    os.makedirs(barcode_output_dir, exist_ok=True)
    output_demuxlet = f"{barcode_output_dir}/demuxlet_output.tsv"

    demuxlet_command = (
        f"singularity exec {popscle_sif} popscle demuxlet "
        f"--plp {plp_file} --field GT --vcf {vcf_file} "
        f"--group-list {full_barcode_list} --out {output_demuxlet}"
    )
    print(f"Running demuxlet command: {demuxlet_command}")

    result = subprocess.run(demuxlet_command, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Error running demuxlet for barcode {barcode}:\n{result.stderr}")
    else:
        print(f"Demuxlet completed successfully for barcode {barcode}")

if __name__ == "__main__":
    with Pool(processes=3) as pool:
        pool.map(run_demuxlet, barcodes)
