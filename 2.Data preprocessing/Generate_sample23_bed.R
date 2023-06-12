library(Signac)
library(Seurat)
library(tidyverse)
library(Matrix)
library(glue)
source('preprocess_functions.R')
sample23_ATAC <- readRDS('../../data/scATAC-seq/Sample23/sample23_scATAC-seq.Rds')

sparse_mtx <- sample23_ATAC@assays$peaks@counts
rownames(sparse_mtx) <- paste0('chr',rownames(sparse_mtx))

##subset columns, if your peak is less than 1M, your does not need to use multithreads

pseudobulk <- generate_pseudobulk(sparse_mtx,group_by = sample23_ATAC$cell_type)

cell_gr <- separate(as.data.frame(rownames(pseudobulk)),col = everything(),sep = '-',into = c('seqnames','start','end'))
cell_gr_list <- map(colnames(pseudobulk),get_cell_gr,seurat_Obj = sample23_ATAC)
names(cell_gr_list) <- colnames(pseudobulk)

sample <- 'sample23'

for (i in 1:length(cell_gr_list)) {
  if (!dir.exists(glue('../../data/bed/{sample}'))) {
    dir.create(glue('../../data/bed/{sample}'))
  }
  filenames <- paste0(glue('../../data/bed/{sample}/'),names(cell_gr_list)[i],'.bed.gz')
  rtracklayer::export.bed(object = cell_gr_list[[i]],con = filenames)
  cat(names(cell_gr_list)[i],'\n')
}
