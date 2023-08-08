library(Seurat)
library(Signac)
library(sceasy)
library(tidyverse)
library(reticulate)
use_condaenv('scanpy')
loompy <- reticulate::import('loompy')
source('../2.Data preprocessing/preprocess_functions.R')
#sample1----
sample1_ATAC <- readRDS('../../data/scATAC-seq/sample1/Rds/sample1_scATAC-seq_80k_processed.Rds')
sceasy::convertFormat(sample1_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/sample1/Rds/sample1_scATAC-seq_80k_processed.h5ad')
#sample4----
sample4_mtx <- readRDS('../../data/scATAC-seq/Sample4/Sample4_peak_matrix.Rds')
sparse_mtx <- sample4_mtx@assays@data$PeakMatrix
rownames(sparse_mtx) <- paste(as.data.frame(sample4_mtx@rowRanges)[[1]],as.data.frame(sample4_mtx@rowRanges)[[2]],as.data.frame(sample4_mtx@rowRanges)[[3]],sep = '-')
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sparse_mtx),subset_peaks)),]
sample4_mtx$cell_type <- sample4_mtx$Celltype
chrom_assay <- CreateChromatinAssay(
  counts = sparse_mtx,
  sep = c("-", "-"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)
sample4_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)
identical(colnames(sparse_mtx),rownames(sample4_mtx@colData))
sample4_ATAC <- AddMetaData(sample4_ATAC,metadata = as.data.frame(sample4_mtx@colData))
saveRDS(sample4_ATAC,file = '../../data/scATAC-seq/Sample4/sample4_scATAC-seq_30k_processed.Rds')
sceasy::convertFormat(sample4_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/sample4/sample4_scATAC-seq_30k_processed.h5ad')
#sample7----
sample7_mtx<- readRDS('../../data/scATAC-seq/Sample7/Sample7_peak_matrix.Rds')
sparse_mtx <- sample7_mtx@assays@data$PeakMatrix
sample7_mtx$cell_type <- sample7_mtx$Sample
rownames(sparse_mtx) <- paste(as.data.frame(sample7_mtx@rowRanges)[[1]],as.data.frame(sample7_mtx@rowRanges)[[2]],as.data.frame(sample7_mtx@rowRanges)[[3]],sep = '-')
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sparse_mtx),subset_peaks)),]
chrom_assay <- CreateChromatinAssay(
  counts = sparse_mtx,
  sep = c("-", "-"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)
sample7_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = as.data.frame(sample7_mtx@colData)
)

saveRDS(sample7_ATAC,file = '../../data/scATAC-seq/Sample7/sample7_scATAC-seq_50k_processed.Rds')

sceasy::convertFormat(sample7_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/sample7/sample7_scATAC-seq_50k_processed.h5ad')
