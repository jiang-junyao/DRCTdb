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
  if (is.na(str_extract(x,'.*(?=\\s\\d)'))) {
    res <- x
  }else{
    res <- str_extract(x,'.*(?=\\s\\d)')
  }
}))

merged_cell_type <- unique(cell_type$new_cell_type)


reduced_cell_type <- map(seq_along(merged_cell_type),function(i){
  concensus_gr <- GenomicRanges::reduce(reduce(cell_gr_list[str_which(names(cell_gr_list),fixed(merged_cell_type[i]))],c))
  return(concensus_gr)
}) %>% setNames(merged_cell_type)


names(reduced_cell_type) <- str_replace(names(reduced_cell_type),'\\+','_')

for (i in 1:length(reduced_cell_type)) {
  filenames <- paste0('../../data/bed/sample3/',names(reduced_cell_type)[i],'.bed.gz')
  rtracklayer::export.bed(object = reduced_cell_type[[i]],con = filenames)
  cat(names(reduced_cell_type)[i],'\n')
}


