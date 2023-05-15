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



cell_type <- tibble(
  raw_cell_type = names(cell_gr_list),
  new_cell_type = map_chr(names(cell_gr_list),function(x){
    res <- str_replace_all(x,' ', '_')
    return(res)
}))

names(cell_gr_list) <- cell_type$new_cell_type



  
sample3 <- cell_gr_list[-str_which(names(cell_gr_list),'Fetal')]



for (i in 1:length(sample3)) {
  filenames <- paste0('../../data/bed/sample3/',names(sample3)[i],'.bed.gz')
  rtracklayer::export.bed(object = sample3[[i]],con = filenames)
  cat(names(sample3)[i],'\n')
}

sample5 <- cell_gr_list[str_which(names(cell_gr_list),'Fetal')]

for (i in 1:length(sample5)) {
  filenames <- paste0('../../data/bed/sample5/',names(sample5)[i],'.bed.gz')
  rtracklayer::export.bed(object = sample5[[i]],con = filenames)
  cat(names(sample5)[i],'\n')
}

