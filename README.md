# DRCTdb (temporary)
Disease related cell type analysis to decode cell type effect and underlying regulatory mechanisms across human disease


1. Formatting data from different sources

2. Preprocessing data to get following ouputs: (I) scRNA-seq: seurat object with cell annotations in metadata; (II) scATAC-seq: seurat or ArchR object with cell annotations in metadata; (III) eqtl: SNPs region and related gene expression.

3. Integrating scATAC-seq & GWAS to get disease related cell types.

4. Integrating scRNA-seq & eQTL to identify disease related cell type marker genes.

5. Downstream analysis: (I) Use disease related marker genes to perform co-expression analysis; (II) TF-target analysis in marker regions; (III) Integrate information from I and II to predict underlying regulation mechanisms in disease related cell types. (IV) identify which SNPs open chromatin to cause disease and which SNPs close chromatin to cause disease
