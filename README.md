# DRCTdb (temporary)
Disease related cell type analysis to decode cell type effect and underlying regulatory mechanisms across human disease


1. Formatting data from different sources

2. Preprocessing data to get following ouputs: 

- scRNA-seq: seurat object with cell annotations in metadata
- scATAC-seq: seurat or ArchR object with cell annotations in metadata
- eQTL: SNPs region and related gene expression.

3. Calculate the cell specific genomic intervals 

- scRNAseq:
Top cell specific expressed gene(10%) based on the t statitic model (Described detail in [Finucane et al. 2018 Nat Genet.](https://www.nature.com/articles/s41588-018-0081-4) ), which used to determine the genetic window (100kb) as input.

- scATAC-seq:
Cell specific open chromtin regions were uesed as input for LDSC

    If fine mapped cell specific cCREs available, used cCREs replace to Peaks.

4. Using [LDSC](https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses) to identify disease related cell type marker genes.

5. Downstream analysis: 
- (I) Use disease related marker genes to perform co-expression analysis; 
- (II) TF-target analysis in marker regions; 
- (III) Integrate information from I and II to predict underlying regulation mechanisms in disease related cell types.
- (IV) Identify which SNPs open chromatin to cause disease and which SNPs close chromatin to cause disease
