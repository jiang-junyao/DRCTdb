library(Signac)
library(Seurat)
library(tidyverse)
library(Matrix)
library(glue)
source('preprocess_functions.R')
sample1_ATAC <- readRDS('../../data/scATAC-seq/sample16/sample16_Bone_marrow_ATAC_Healthy_35k_processed.Rds')

sparse_mtx <- sample1_ATAC@assays$peaks@counts
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sparse_mtx),subset_peaks)),]

pseudobulk <- generate_pseudobulk(sparse_mtx,group_by = sample1_ATAC$cell_type)
cell_gr <- tidyr::separate(as.data.frame(rownames(pseudobulk)),col = everything(),sep = '-',into = c('seqnames','start','end'))

cell_gr_list <- purrr::map(pseudobulk,function(x){
  threshold <- str_extract(colnames(as.data.frame(catable(x)))[5],'\\d+') %>% as.numeric()
  gr <- cell_gr[which(as.vector(x) >= threshold),] %>% GenomicRanges::makeGRangesFromDataFrame()
  return(gr)
})

sample <- 'sample16'

for (i in 1:length(cell_gr_list)) {
  if (!dir.exists(glue('../../data/bed/{sample}'))) {
    dir.create(glue('../../data/bed/{sample}'))
  }
  filenames <- paste0(glue('../../data/bed/{sample}/'),names(cell_gr_list)[i],'.bed.gz')
  rtracklayer::export.bed(object = cell_gr_list[[i]],con = filenames)
  cat(names(cell_gr_list)[i],'\n')
}

