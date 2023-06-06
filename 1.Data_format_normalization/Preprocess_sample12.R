library(Seurat)
library(Signac)
load('../../data/scATAC-seq/Sample12/SNARE-Seq2/GSE183273_Kidney_Healthy-Injury_Cell_Atlas_SNARE2-RNA-AC_Seurat_03282022.rda')
metadata <-  KSAC@meta.data
sample12_RNA <- CreateSeuratObject(KSAC@assays$RNA@counts,meta.data = KSAC@meta.data , project = 'sample12_RNA')

sample12_RNA$cell_type <- sample12_RNA$subclass.full
saveRDS(sample12_RNA,'../../data/scRNA-seq/Sample12/Sample12_RNA-seq.Rds')

chrom_assay <- CreateChromatinAssay(
  counts = KSAC@assays$ATAC$counts,
  sep = c("-", "-"),
  genome = 'hg38',
  min.cells = 1,
  min.features = 1
)
sample12_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = KSAC@meta.data 
)
sample12_ATAC$cell_type <- sample12_ATAC$subclass.full
saveRDS(sample12_ATAC,'../../data/scATAC-seq/Sample12/Sample12_ATAC-seq.Rds')
