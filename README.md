# DRCTdb (temporary)
Disease related cell type analysis to decode cell type effect and underlying regulatory mechanisms across human disease


1. Formatting data from different sources

2. Preprocessing data to get following ouputs: (I) scRNA-seq: seurat object with cell annotations in metadata; (II) scATAC-seq: seurat or ArchR object with cell annotations in metadata; (III) eqtl: SNPs region and related gene expression.

3. Integrating scRNA-seq and eqtl data to get disease related cell types and marker genes

4. Integrating scATAC-seq data to get disease related cell types and marker regions

5. Machine learning model to calculate weight of cell types related to disease.

6. Downstream analysis: (I) Use disease related marker genes to perform co-expression analysis; (II) TF-target analysis in marker regions; (III) Integrate information from I and II to predict underlying regulation mechanisms in disease related cell types.
