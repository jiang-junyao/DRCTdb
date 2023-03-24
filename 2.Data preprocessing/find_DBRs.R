### Preprocess scATAC-seq to get DBRs
find_DBRs <- function(obj,test_method = 'wilcox'){
  library(Signac)
  library(GenomicRanges)
  obj <- RunTFIDF(obj)
  DBRs <- FindAllMarkers(obj,test.use = test_method)
  DBRs_gr <- GRanges(rownames(DBRs))
  DBRs_gr$cluster <- DBRs$cluster
  DBRs_gr$p_adj <- DBRs$p_val_adj
  DBRs_gr$log2fc <- DBRs$avg_log2FC
  return(DBRs_gr)
}




