import subprocess
import os
import sys

index = int(sys.argv[1])
cell_barcodes_path = f"barcodes/cell_barcodes_{index}.txt"
output_directory = "bams_Flu.1_new"
os.makedirs(output_directory, exist_ok=True)
input_bam = "../alignment/output_of_alignment_YPD2_with_clean_gtf/run_count/outs/possorted_genome_bam.bam"
tmpdir = os.environ.get('TMPDIR', '/tmp')

with open(cell_barcodes_path, 'r') as f:
    barcodes = [line.strip() for line in f if line.strip()]

all_barcodes_file = os.path.join(tmpdir, f"barcodes_{index}.txt")
with open(all_barcodes_file, 'w') as f:
    f.write('\n'.join(barcodes))

subset_bam = os.path.join(tmpdir, f"subset_{index}.bam")
subprocess.run(f"/projects/rps/cgsb/gencore/mk5636/yeast_strain_analysis/subset-bam_linux --bam {input_bam} --cell-barcodes {all_barcodes_file} --out-bam {subset_bam} --cores 2", shell=True, check=True)

sorted_bam = os.path.join(tmpdir, f"subset_{index}_sorted.bam")
subprocess.run(f"samtools sort -o {sorted_bam} {subset_bam}", shell=True, check=True)
subprocess.run(f"samtools index {sorted_bam}", shell=True, check=True)
os.remove(subset_bam)

for barcode in barcodes:
    output_bam = os.path.join(output_directory, f"{barcode}.bam")
    subprocess.run(f"samtools view -b -d CB:{barcode} {sorted_bam} > {output_bam}", shell=True, check=True)

os.remove(sorted_bam)
os.remove(f"{sorted_bam}.bai")
os.remove(all_barcodes_file)
