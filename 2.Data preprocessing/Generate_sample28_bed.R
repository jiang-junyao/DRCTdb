library(Signac)
library(Seurat)
library(tidyverse)
library(Matrix)
library(glue)
source('preprocess_functions.R')
sample28_ATAC <- readRDS('../../data/scATAC-seq/Sample28/sample28_scATAC-seq_29k_processed.Rds')
sample28_ATAC$cell_type <- str_replace_all(sample28_ATAC$cell_type,' ','_')
sample28_ATAC$cell_type <- str_replace_all(sample28_ATAC$cell_type,'/','_')
sparse_mtx <- sample28_ATAC@assays$peaks@counts
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sample28_ATAC),subset_peaks)),]

pseudobulk <- generate_pseudobulk(sparse_mtx,group_by = sample28_ATAC$cell_type)
cell_gr <- tidyr::separate(as.data.frame(rownames(pseudobulk)),col = everything(),sep = '-',into = c('seqnames','start','end'))
cell_gr_list <- map(colnames(pseudobulk),get_cell_gr,seurat_Obj = sample28_ATAC)
names(cell_gr_list) <- colnames(pseudobulk)
sample <- 'sample28'
for (i in 1:length(cell_gr_list)) {
    if (!dir.exists(glue('../../data/bed/{sample}'))) {
        dir.create(glue('../../data/bed/{sample}'))
    }
    filenames <- paste0(glue('../../data/bed/{sample}/'),names(cell_gr_list)[i],'.bed.gz')
    rtracklayer::export.bed(object = cell_gr_list[[i]],con = filenames)
    cat(names(cell_gr_list)[i],'\n')
}
