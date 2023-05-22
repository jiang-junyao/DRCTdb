library(Seurat)

sample8_ATAC <- readRDS('../../data/scATAC-seq/Sample8/local (1).rds')
atac_metadata <- sample8_ATAC@meta.data
data.table::fwrite(atac_metadata,'../../data/scATAC-seq/Sample8/sample8_ATAC_metadata.txt',sep = '\t')

sample8_RNA <- readRDS('../../data/scRNA-seq/Sample8/local.rds')
RNA_metadata <- sample8_RNA@meta.data
data.table::fwrite(RNA_metadata,'../../data/scATAC-seq/Sample8/sample8_RNA_metadata.txt',sep = '\t')