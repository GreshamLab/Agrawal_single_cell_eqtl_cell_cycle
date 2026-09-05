import scanpy as sc

# Load the ALL-GENES saved adata (not the program-select one)
adata = sc.read_h5ad(
    '/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all'
    '/YPD.1/continuous_time_trajectory/data/adata_cc_time_all_genes_spellman_annotated.h5ad'
)

print(adata)
print(adata.obs.columns.tolist())      # confirm 'CC' and 'cc_time' are here

# Export obs to CSV — new analysis has ONE time column: cc_time
adata.obs[['CC', 'cc_time']].to_csv(
    '/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all'
    '/YPD.1/continuous_time_trajectory/data/cc_time_all_genes_spellman_annotated.csv'
)
print("CSV saved!")
print(adata.obs[['CC', 'cc_time']].head())
