library(Signac)
library(Seurat)
library(tidyverse)
library(Matrix)
library(glue)
source('preprocess_functions.R')
sample1_ATAC <- readRDS('../../data/scATAC-seq/sample1/Rds/sample1_scATAC-seq_80k_processed.Rds')

sparse_mtx <- sample1_ATAC@assays$peaks@counts
sparse_mtx <- sparse_mtx[str_starts(rownames(sparse_mtx),'chr'),]

##subset columns, if your peak is less than 1M, your does not need to use multithreads

pseudobulk <- generate_pseudobulk(sparse_mtx,group_by = sample1_ATAC$cell_type)

# pseudobulk <- generate_pseudobulk(sparse_mtx,group_by = na.omit(sample1_ATAC$cell_type))
# identical(rownames(pseudobulk),rownames(sparse_mtx))
# system.time(generate_pseudobulk(sparse_mtx,group_by = na.omit(sample1_ATAC$cell_type)))
# 
# pseudobulk2 <- generate_pseudobulk(sparse_mtx,group_by = na.omit(sample1_ATAC$cell_type),ncore = 10)
# identical(rownames(pseudobulk),rownames(sparse_mtx))
# identical(rownames(pseudobulk),rownames(pseudobulk2))
# system.time(generate_pseudobulk(sparse_mtx,group_by = na.omit(sample1_ATAC$cell_type),ncore = 10))

cell_gr <- separate(as.data.frame(rownames(pseudobulk)),col = everything(),sep = '-',into = c('seqnames','start','end'))




cell_gr_list <- purrr::map(pseudobulk,function(x){
  threshold <- str_extract(colnames(as.data.frame(catable(x)))[4],'\\d+') %>% as.numeric()
  return(gr)
})


cell_gr_list <- map(colnames(pseudobulk),get_cell_gr,seurat_Obj = sample1_ATAC)
names(cell_gr_list) <- colnames(pseudobulk)


sample <- 'sample1'

for (i in 1:length(cell_gr_list)) {
  if (!dir.exists(glue('../../data/bed/{sample}'))) {
    dir.create(glue('../../data/bed/{sample}'))
  }
  filenames <- paste0(glue('../../data/bed/{sample}'),names(cell_gr_list)[i],'.bed.gz')
  rtracklayer::export.bed(object = cell_gr_list[[i]],con = filenames)
  cat(names(cell_gr_list)[i],'\n')
}

