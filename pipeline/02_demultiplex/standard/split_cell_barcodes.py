def split_file(filename, lines_per_file, output_prefix):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    total_files = len(lines) // lines_per_file
    if len(lines) % lines_per_file != 0:
        total_files += 1
    for i in range(total_files):
        with open(f"{output_prefix}_{i}.txt", 'w') as file:
            file.writelines(lines[i*lines_per_file:(i+1)*lines_per_file])

# Call the function to split 'cell_barcodes.txt' into files with prefix 'cell_barcodes'
split_file('../gene_expression_analysis/cell_barcodes_Q4.txt', 536, 'barcodes_Q4/cell_barcodes')

