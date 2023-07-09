library(Seurat)
sample8_scRNA <- readRDS('../../data/scRNA-seq/Sample8/local.rds')
saveRDS(sample8_scRNA,file = '../../data/scRNA-seq/Sample8/sample8_scRNA_processed_68k.Rds')


metadata <- data.table::fread('../../data/scATAC-seq/Sample8/sample8_ATAC_metadata.txt')
