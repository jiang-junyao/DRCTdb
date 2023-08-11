### Transfer data to seurat object
library(Seurat)
library(Signac)
library(tidyverse)
source('../2.Data preprocessing/preprocess_functions.R')
sample23_mtx <- readRDS('../../data/scATAC-seq/sample23/islet_reproducible_peak_matrix.Rds')

sparse_mtx <- sample23_mtx@assays@data$PeakMatrix
rownames(sparse_mtx) <- paste(as.data.frame(sample23_mtx@rowRanges)[[1]],as.data.frame(sample23_mtx@rowRanges)[[2]],as.data.frame(sample23_mtx@rowRanges)[[3]],sep = '-')
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sparse_mtx),subset_peaks)),]
chrom_assay <- CreateChromatinAssay(
  counts = sparse_mtx,
  sep = c("-", "-"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)
sample23_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = as.data.frame(sample23_mtx@colData)
)
saveRDS(sample23_ATAC,file = '../../data/scATAC-seq/sample23/sample23_scATAC-seq_95k_processed.Rds')



##

