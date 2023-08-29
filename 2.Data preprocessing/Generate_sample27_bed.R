library(Signac)
library(Seurat)
library(tidyverse)
library(Matrix)
library(glue)
source('preprocess_functions.R')
sample27_ATAC <- readRDS('../../data/scATAC-seq/Sample27/Sample27_scATAC-seq_133k_processed.Rds')
sample27_ATAC$cell_type <- str_replace_all(sample27_ATAC$cell_type,' ','_')
sample27_ATAC$cell_type <- str_replace_all(sample27_ATAC$cell_type,'/','_')

sparse_mtx <- sample27_ATAC@assays$peaks@counts
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sample27_ATAC),subset_peaks)),]

pseudobulk <- generate_pseudobulk(sparse_mtx,group_by = sample27_ATAC$cell_type)
cell_gr <- tidyr::separate(as.data.frame(rownames(pseudobulk)),col = everything(),sep = '-',into = c('seqnames','start','end'))
cell_gr_list <- map(colnames(pseudobulk),get_cell_gr,seurat_Obj = sample27_ATAC)
names(cell_gr_list) <- colnames(pseudobulk)
sample <- 'sample27'

for (i in 1:length(cell_gr_list)) {
    if (!dir.exists(glue('../../data/bed/{sample}'))) {
        dir.create(glue('../../data/bed/{sample}'))
    }
    filenames <- paste0(glue('../../data/bed/{sample}/'),names(cell_gr_list)[i],'.bed.gz')
    rtracklayer::export.bed(object = cell_gr_list[[i]],con = filenames)
    cat(names(cell_gr_list)[i],'\n')
}