# DRCTdb
Disease related cell type analysis to decode cell type effect and underlying regulatory mechanisms across human disease

1. Formatting data from different sources

2. Preprocessing data to get the following ouputs: 

- scRNA-seq: seurat object with cell annotations in metadata
- scATAC-seq: seurat or ArchR object with cell annotations in metadata

3. Calculate the cell specific genomic intervals by scATAC-seq

    3.1 Call Cell specific open chromatin regions were used as input for LDSC

    3.2 Using [LDSC](https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses) to identify disease related cell type marker genes.

    3.3 Downstream analysis: 
   - (I) Overlap disease related SNPs and cell-type specific variable CREs. In disease related samples, we suggest these SNPs open chromatin to cause disease. In normal sample, we suggest these SNPSs close chromatin to cause disease.
   - (II) (Optional) Identify cell type specific variable genes that are related SNPs by overlapping the gene's TSS (Transcription Start Site) region with the SNPs.
   - (III) Perform a Transcription Factor (TF) targeting analysis in the CREs identified in part I (significant motif related tf as tf, region annotated gene as target). This will help us further understand how the TFs might play a role in regulating specific genes in respective cell types.
   - (IV) When scATAC-seq data is available, integrating Pearson correlation and TF-target result of part IIII to reconstruct cell type-specific gene regulatory network. When scATAC-seq data is not available, using [SCENIC](https://scenicplus.readthedocs.io/en/latest/index.html) to reconstruct cell type-specific diseased related gene regulatory network.


We offer a step-by-step tutorials to preprocess single cell multiomics data in **reproduce_case folder**

DRCTDB is a fully open-source database, where users can download the entire project and run locally at **db-core folder**