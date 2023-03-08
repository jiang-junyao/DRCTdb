### Transfer data to seurat object
library(Seurat)
library(GenomicRanges)
library(tidyverse)
merged_matrix <- ReadMtx(features = '../../data/bed/sample3/rawdata/cCREs.bed.gz',
                         cells = '../../data/bed/sample3/rawdata/celltypes.txt.gz',
                         mtx = '../../data/bed/sample3/rawdata/matrix.tsv.gz',
                         feature.column = 1)
cres <- data.table::fread('../../data/bed/sample3/rawdata/cCREs.bed.gz')
colnames(cres) <- c('seqnames','start','end')

merged_matrix_df <- as.data.frame(merged_matrix)

cell_gr_list <- purrr::map(merged_matrix_df,function(x){
  gr <- cres[which(x == 1)] %>% makeGRangesFromDataFrame()
  return(gr)
})

for (i in 1:length(cell_gr_list)) {
  filenames <- paste0('../../data/bed/sample3/',names(cell_gr_list)[i],'.bed.gz')
  rtracklayer::export.bed(object = cell_gr_list[[i]],con = filenames)
  cat(names(cell_gr_list)[i],'\n')
}


