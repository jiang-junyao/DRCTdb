library(Signac)
library(Seurat)
library(tidyverse)
library(Matrix)
library(glue)
source('preprocess_functions.R')
sample8_ATAC <- readRDS('../../data/scATAC-seq/Sample8/sample8_peak_matrix.Rds')
sparse_mtx <- sample8_ATAC@assays@data$PeakMatrix
rownames(sparse_mtx) <- paste(as.data.frame(sample8_ATAC@rowRanges)[[1]],as.data.frame(sample8_ATAC@rowRanges)[[2]],as.data.frame(sample8_ATAC@rowRanges)[[3]],sep = '-')


pseudobulk <- generate_pseudobulk(sparse_mtx,group_by = sample8_ATAC$cell_type)

cell_gr <- tidyr::separate(as.data.frame(rownames(pseudobulk)),col = everything(),sep = '-',into = c('seqnames','start','end'))
cell_gr_list <- map(colnames(pseudobulk),get_cell_gr,seurat_Obj = sample8_ATAC,.progress = T)
names(cell_gr_list) <- colnames(pseudobulk)

sample <- 'sample8'


for (i in 1:length(cell_gr_list)) {
  if (!dir.exists(glue('../../data/bed/{sample}'))) {
    dir.create(glue('../../data/bed/{sample}'))
  }
  filenames <- paste0(glue('../../data/bed/{sample}/'),names(cell_gr_list)[i],'.bed.gz')
  rtracklayer::export.bed(object = cell_gr_list[[i]],con = filenames)
  cat(names(cell_gr_list)[i],'\n')
}
