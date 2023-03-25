### Preprocess scATAC-seq to get DBRs
find_DBRs <- function(obj,test_method = 'wilcox'){
  library(Signac)
  library(GenomicRanges)
  obj <- RunTFIDF(obj)
  DBRs <- FindAllMarkers(obj,test.use = test_method)
  str_peak <- length(strsplit(DBRs$gene[1],'-')[[1]])
  if (str_peak == 2) {
    DBRs_gr <- GRanges(rownames(DBRs))
  }else if(str_peak == 3){
    peak_name = t(as.data.frame(strsplit(DBRs$gene,'-')))
    peak_name = paste0(peak_name[,1],':',peak_name[,2],'-',peak_name[,3])
    DBRs_gr <- GRanges(peak_name)
  }
  
  DBRs_gr$cluster <- DBRs$cluster
  DBRs_gr$p_adj <- DBRs$p_val_adj
  DBRs_gr$log2fc <- DBRs$avg_log2FC
  return(DBRs_gr)
}


