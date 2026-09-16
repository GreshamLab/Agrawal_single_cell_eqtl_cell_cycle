# Single-Cell eQTL Mapping and Cell Cycle Analysis in *S. cerevisiae*
Analysis code and results for a study mapping expression quantitative trait loci (eQTLs) from single-cell RNA-seq data in a genetically diverse panel of wild *Saccharomyces cerevisiae* strains, across growth conditions and cell cycle states.
**Author:** Akriti Agrawal — [Gresham Lab](https://greshamlab.bio.nyu.edu/), New York University

## Overview
Ninety-six strains from the 1011 *S. cerevisiae* genomes collection ([Peter et al. 2018](https://doi.org/10.1038/s41586-018-0030-5)) were pooled and profiled by droplet-based single-cell RNA-seq, with strain identity assigned per cell using genotype demultiplexing (`demuxlet`). 
Pseudobulk expression per strain was used to map cis- and trans-eQTLs (SNPs, indels and structural variants) against strain genotypes, in standard growth medium (YPD) and under two additional conditions (fluconazole treatment and synthetic complete 2% glucose). Cell cycle phase and continuous cell cycle time were also assigned to single cells to examine how eQTL effects vary across the cell cycle

## Repository structure
The repository is organized by manuscript figure. Each figure directory generally contains a `codes/` (or `code/`) subdirectory with the R Markdown / Python / shell scripts used to produce the panel, and a `results/` subdirectory with the resulting plots, tables, and intermediate data files.

```
Figure.1/   Strain panel, demuxlet genotype demultiplexing, and scRNA-seq reproducibility
Figure.2/   Strain- and clade-level differential gene expression
Figure.3/   Genotype preparation, pseudobulk expression, and cis/trans-eQTL mapping in YPD
Figure.4/   eQTL mapping and comparison across conditions (YPD, fluconazole, SC2Glu), plus bulk-vs-single-cell validation
Figure.5/   Cell cycle phase and cell cycle time assignment for single cells, and cell-cycle-resolved phenotype analysis
```
### Figure 1 — Strain panel and scRNA-seq quality control
- **Panel A** — phylogenetic tree of the 1011-strain collection with the 96 sequenced strains highlighted (`Fig.1-Panel-A-strain-selection-tree`).
- **Panel B** — `demuxlet` genotype demultiplexing and validation of strain assignment accuracy against a spiked-in reference strain (S288C) (`Fig.1-Panel-B-demuxlet-accurately-detects-S288C`).
- **Panels C/D** — reproducibility of scRNA-seq expression and strain abundance across replicate runs (`Fig.1-Panel-C-scRNAseq-reproducibility_agreement`).
- Additional experiment summary tables and supplementary QC figures.

### Figure 2 — Differential gene expression by strain and clade
Pseudobulk differential expression between strains and clades, GO over-representation and gene set enrichment analysis, and regulon-level expression (e.g. GAL, NCR, stress response genes), with per-condition supplementary analyses.

### Figure 3 — eQTL mapping in YPD
- Genotype preparation (custom yeast annotations, VCF processing, PLINK LD clumping) for SNPs, indels, structural variants, and CNVs.
- Pseudobulk expression preparation and linear mixed-model eQTL mapping (`fastlmm`/pseudobulk pipeline).
- Cis/trans classification and effect-size analysis of eQTLs by variant type.
- Raw single-cell validation of top eQTL hits (violin/scatter plots per variant).
- Reproducibility of eQTL mapping between biological replicates and conditions.

### Figure 4 — eQTL mapping across conditions
Repeats the Figure 3 eQTL mapping pipeline for fluconazole and SC2Glu conditions, compares eQTL sharing and effect-size concordance across conditions, and validates single-cell (SC2Glu) results against bulk RNA-seq.

### Figure 5 — Cell cycle assignment and cell-cycle-resolved analysis
- Assignment of discrete cell cycle phase to single cells from normalized expression data.
- Assignment of continuous cell cycle time (pseudotime) using cell-cycle marker/harmonic gene sets, including an `inferelator`-based pipeline (Python/bash scripts under `part1-bash-code-to-assign-cell-cycle-time`).
- Identification of genes/strains whose cell-cycle amplitude, timing, or structure vary across the strain panel.

## Requirements

Analyses are implemented primarily in R (R Markdown notebooks), with some Python and shell scripts for cell cycle time inference. Key dependencies include:

- **R**: `Seurat`, `tidyverse` (`dplyr`, `tidyr`, `readr`), `ggplot2`,`ggtree`, `viridis`, `org.Sc.sgd.db`, `Rtsne`
- **Genotyping / eQTL mapping**: `demuxlet`, `PLINK`, `fastlmm`
- **Python**: standard scientific stack (`numpy`, `pandas`, `scanpy`/`h5ad` tooling) plus [`inferelator`](https://github.com/flatironinstitute/inferelator) for cell cycle time inference

Each `.Rmd` file lists its specific library dependencies in its setup chunk. Raw sequencing data and full genotype VCFs are not included in this repository; several notebooks reference external paths to the Gresham Lab compute environment.

## Usage

Each figure's `codes/` scripts are intended to be run in the order implied by their filenames (e.g. `Step1`, `Step2`, ... or `Panel A`, `Panel B`, ...) within that figure's directory, reading inputs from and writing outputs to the
adjacent `results/` (and, where present, `data/`) directories. There is no single top-level pipeline script; each figure is self-contained.

## Citation

If you use this code, please cite the associated manuscript (citation to be added upon publication).
