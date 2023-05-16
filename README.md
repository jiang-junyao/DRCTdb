# DRCTdb (temporary)
Disease related cell type analysis to decode cell type effect and underlying regulatory mechanisms across human disease


1. Formatting data from different sources

2. Preprocessing data to get following ouputs: 

- scRNA-seq: seurat object with cell annotations in metadata
- scATAC-seq: seurat or ArchR object with cell annotations in metadata
- eQTL: SNPs region and related gene expression.

3. Calculate the cell specific genomic intervals 

- scRNAseq:
Top cell specific expressed gene(10%) based on the t statitic model (Described detail in [Finucane et al. 2018 Nat Genet.](https://www.nature.com/articles/s41588-018-0081-4)), which used to determine the genetic window (100kb) as input.

- scATAC-seq:
Cell specific open chromtin regions were uesed as input for LDSC

    If fine mapped cell specific cCREs available, used cCREs replace to Peaks.

4. Using [LDSC](https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses) to identify disease related cell type marker genes.

5. Downstream analysis: 
- (I) Overlap disease related SNPs and cell-type specific variable CREs. In disease related samples, we suggest these SNPs open chromatin to cause disease. In normal sample, we suggest these SNPSs close chromatin to cause disease.
- (II) (Optional) Identify cell type specific variable genes that are related SNPs by overlapping the gene's TSS (Transcription Start Site) region with the SNPs.
- (III) Perform a Transcription Factor (TF) targeting analysis in the CREs identified in part I (significant motif related tf as tf, region annotated gene as target). This will help us further understand how the TFs might play a role in regulating specific genes in respective cell types.
- (IV) When scATAC-seq data is available, integrating Perason correlation and TF-target result of part IIII to reconstruct cell type-specific gene regulatory network. When scATAC-seq data is not available, using [SCENIC](https://scenicplus.readthedocs.io/en/latest/index.html) to reconstruct cell type-specific gene regulatory network.

