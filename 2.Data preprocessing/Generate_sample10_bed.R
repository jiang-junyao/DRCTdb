library(Signac)
library(Seurat)
library(tidyverse)
library(Matrix)
library(glue)
source('preprocess_functions.R')
sample10_ATAC <- readRDS('../../data/scATAC-seq/Sample10/sample10_peak_matrix.Rds')
sample10_ATAC$cell_type <- str_replace(sample10_ATAC$cell_type,'\\+','_high')
sample10_ATAC$cell_type <- str_replace_all(sample10_ATAC$cell_type,'/','_')
sample10_ATAC$cell_type <- str_replace_all(sample10_ATAC$cell_type,' ','_')
sparse_mtx <- sample10_ATAC@assays@data$PeakMatrix
rownames(sparse_mtx) <- paste(as.data.frame(sample10_ATAC@rowRanges)[[1]],as.data.frame(sample10_ATAC@rowRanges)[[2]],as.data.frame(sample10_ATAC@rowRanges)[[3]],sep = '-')
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sparse_mtx),subset_peaks)),]

pseudobulk <- generate_pseudobulk(sparse_mtx,group_by = sample10_ATAC$cell_type)
cell_gr <- tidyr::separate(as.data.frame(rownames(pseudobulk)),col = everything(),sep = '-',into = c('seqnames','start','end'))
cell_gr_list <- map(colnames(pseudobulk),get_cell_gr,seurat_Obj = sample10_ATAC)
names(cell_gr_list) <- colnames(pseudobulk)

sample <- 'sample10'

for (i in 1:length(cell_gr_list)) {
  if (!dir.exists(glue('../../data/bed/{sample}'))) {
    dir.create(glue('../../data/bed/{sample}'))
  }
  filenames <- paste0(glue('../../data/bed/{sample}/'),names(cell_gr_list)[i],'.bed.gz')
  rtracklayer::export.bed(object = cell_gr_list[[i]],con = filenames)
  cat(names(cell_gr_list)[i],'\n')
}
