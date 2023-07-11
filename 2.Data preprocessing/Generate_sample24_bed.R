library(Signac)
library(Seurat)
library(tidyverse)
library(Matrix)
library(glue)
source('preprocess_functions.R')
###
sample24_ATAC <- readRDS('../../data/scATAC-seq/Sample24/sample24_scATAC-seq_100k_healthy_processed.Rds')

sparse_mtx <- sample24_ATAC@assays$peaks@counts
sparse_mtx <- sparse_mtx[str_starts(rownames(sparse_mtx),'chr'),]

pseudobulk <- generate_pseudobulk(sparse_mtx,group_by = sample24_ATAC$cell_type)
cell_gr <- separate(as.data.frame(rownames(pseudobulk)),col = everything(),sep = '-',into = c('seqnames','start','end'))

cell_gr_list <- map(colnames(pseudobulk),get_cell_gr,seurat_Obj = sample24_ATAC)
names(cell_gr_list) <- colnames(pseudobulk)

sample <- 'sample24'

for (i in 1:length(cell_gr_list)) {
  if (!dir.exists(glue('../../data/bed/{sample}'))) {
    dir.create(glue('../../data/bed/{sample}'))
  }
  filenames <- paste0(glue('../../data/bed/{sample}/'),names(cell_gr_list)[i],'.bed.gz')
  rtracklayer::export.bed(object = cell_gr_list[[i]],con = filenames)
  cat(names(cell_gr_list)[i],'\n')
}
