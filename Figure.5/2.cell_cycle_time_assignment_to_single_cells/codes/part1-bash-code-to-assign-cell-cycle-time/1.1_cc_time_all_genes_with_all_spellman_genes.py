"""
Cell-cycle pseudotime WITHOUT program_select, on ALL genes.
Standalone analysis: no program discovery, no gene subsetting.
calculate_times is run on the full count matrix.
(The program_select-based analysis is kept separately; this is a new,
independent time assignment on everything.)
"""
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import zscore
from inferelator_velocity.times import calculate_times

ROOT = ('/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all'
        '/YPD.1')

# ── Step 1: Load raw 10x counts ───────────────────────────────────────────────
adata = sc.read_10x_h5(
    f'{ROOT}/alignment/output_of_alignment_YPD1_with_clean_gtf'
    '/run_count/outs/filtered_feature_bc_matrix.h5'
)
adata.var_names_make_unique()
adata.layers['counts'] = adata.X.astype(np.int64)
adata.obs['n_counts'] = np.asarray(adata.layers['counts'].sum(axis=1)).ravel()
print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")

# ── Step 2: Spellman phase membership — match on gene_ids (ORF) ───────────────
cc = pd.read_csv(
    f'{ROOT}/cell_cycle_fastlmm/data/cell_cycle_genes_orf_and_common_names.csv'
)
print(cc.head())
print(f"Spellman genes: {cc['Gene'].nunique()}")

CC_COLS = ['M-G1', 'G1', 'S', 'G2', 'M']

# var_names is 'common name where one exists, else ORF'; gene_ids is ALWAYS ORF.
# Spellman's 'Gene' column is ORFs, so match against gene_ids.
if 'gene_ids' in adata.var.columns:
    orf = adata.var['gene_ids'].astype(str)
else:
    print("WARNING: no 'gene_ids' column — falling back to var_names")
    orf = pd.Series(adata.var_names, index=adata.var_names, dtype=str)

for phase in CC_COLS:
    phase_orfs = cc.loc[cc['Group'] == phase, 'Gene'].dropna().unique()
    adata.var[phase] = orf.isin(phase_orfs).values
    print(f"{phase:5s}: {int(adata.var[phase].sum()):3d} / {len(phase_orfs):3d} matched (on gene_ids/ORF)")

n_matched = int(adata.var[CC_COLS].any(axis=1).sum())
print(f"\nTotal Spellman genes matched: {n_matched} / {cc['Gene'].nunique()}")

multi = adata.var[CC_COLS].sum(axis=1) > 1
print(f"Genes assigned to >1 phase: {int(multi.sum())}")
if multi.sum():
    print(adata.var.loc[multi, CC_COLS])

# ── Step 3: Assign each cell a phase (Spellman-gene z-score argmax) ────────────
def _call_cc(data):
    prop = np.zeros((data.shape[0], len(CC_COLS)), order="F")
    for i, c in enumerate(CC_COLS):
        prop[:, i] = (
            np.asarray(data.layers['counts'][:, data.var[c].values].sum(axis=1)).ravel()
            / data.obs['n_counts'].values
        )
    prop = zscore(prop, ddof=1)
    data.obs['CC'] = np.array(CC_COLS)[np.argmax(prop, axis=1)]
    data.obs['CC'] = data.obs['CC'].astype(
        pd.CategoricalDtype(categories=CC_COLS, ordered=True))

_call_cc(adata)
print("Phase counts:\n", adata.obs['CC'].value_counts())

# ── Step 4: Anchor times (minutes), 88-min wrapped cycle ──────────────────────
CC_TIME_ORDER = {
    'M-G1': ('G1',   7,   22.5),
    'G1':   ('S',    22.5, 39.5),
    'S':    ('G2',   39.5, 56.5),
    'G2':   ('M',    56.5, 77.5),
    'M':    ('M-G1', 77.5, 95),
}

# ── Step 5: ONE calculate_times call on ALL GENES ─────────────────────────────
times, pcs, uns = calculate_times(
    adata.layers['counts'],          # full count matrix — every gene
    adata.obs['CC'].values,
    CC_TIME_ORDER,
    return_components=True,
    verbose=True,
    n_comps=None,                    # MCV/CV picks the number of PCs for the projection
    wrap_time=88.0,
    standardization_method='log_scale',
    mcv_kwargs={'standardization_method': 'log_scale'},
)
adata.obs['cc_time'] = times
adata.obsm['cc_pca'] = pcs
adata.uns['cc_time'] = uns
print("\ncc_time summary:\n", adata.obs['cc_time'].describe())

# ── Step 6: Monotonicity check (expect M-G1 < G1 < S < G2 < M) ────────────────
mean_t = adata.obs.groupby('CC', observed=True)['cc_time'].mean().reindex(CC_COLS)
print("\nMean cc_time per phase:")
print(mean_t)
diffs = np.diff(mean_t.values)
print("consecutive diffs:", np.round(diffs, 1),
      "| monotonic increasing:", bool(np.all(diffs > 0)))

# ── Step 7: Save ──────────────────────────────────────────────────────────────
adata.write_h5ad(f'{ROOT}/continuous_time_trajectory/data/adata_cc_time_all_genes_spellman_annotated.h5ad')
print("\nSaved adata_cc_time_all_genes_spellman_annotated.h5ad")
