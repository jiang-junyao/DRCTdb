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

#Only select peaks fomart as chrxx-xx-xxx
subset_peaks <- function(x){
  if (length(str_extract_all(x,'-')[[1]])== 2) {
    if (str_starts(x,'chr')) {
      return(TRUE)
    }else{
      return(FALSE)
    }
    return(TRUE)
  }else{
    return(FALSE)
  }
}
### Preprocess scATAC-seq to get generate pseudobulk peaks
#input sparse peaks by cell matrix
generate_pseudobulk <- function(mat,group_by,ncore = NULL){
  cell_type  <-  unique(na.omit(group_by))
  row_name <- rownames(mat)
  if (is.null(ncore)) {
    result <- map_dfc(cell_type,function(x,mat){
      col_index <- which(group_by == x)  
      mat <- mat[,col_index]
      mat_rowsum <- as.data.frame(Matrix::rowSums(mat)) 
      colnames(mat_rowsum) <- x
      return(mat_rowsum)
    },mat = mat)
  }else{
    require(furrr)
    plan(multisession, workers = ncore)
    result <- furrr::future_map_dfc(cell_type,function(x,mat){
      col_index <- which(colnames(mat) == x)  
      mat <- mat[,col_index]
      mat_rowsum <- as.data.frame(Matrix::rowSums(mat)) 
      colnames(mat_rowsum) <- x
      return(mat_rowsum)
    },mat = mat)
  }
  return(result)
}


catable <- function (data, categories = c(quantile(data, c(0.01, 0.1, 0.5, 0.9, 0.99), na.rm = TRUE)), cumulative = FALSE, na.rm = TRUE, digits = 3){
  if (!is(data, "numeric")) 
    stop("data should be numeric vector")
  if (!is(categories, "numeric")) 
    stop("categories should be numeric vector")
  ouv <- rep(NA, length(data))
  categories <- sort(categories)
  outmat <- matrix(rep(0, 2 * (length(categories) + 1)), nrow = 2)
  tot <- sum(!is.na(data), na.rm = na.rm)
  outmat[1, 1] <- sum(data <= categories[1], na.rm = na.rm)
  outmat[2, 1] <- (outmat[1, 1]/tot)
  for (i in 1:(length(categories) - 1)) {
    outmat[1, i + 1] <- sum(data > categories[i] & data <= 
                              categories[i + 1], na.rm = na.rm)
    outmat[2, i + 1] <- (outmat[1, i + 1]/tot)
  }
  outmat[1, length(categories) + 1] <- sum(data > categories[length(categories)], 
                                           na.rm = na.rm)
  outmat[2, length(categories) + 1] <- (outmat[1, length(categories) + 
                                                 1]/tot)
  if (cumulative) {
    for (i in 2:(length(categories) + 1)) {
      outmat[1, i] <- (outmat[1, i] + outmat[1, i - 1])
      outmat[2, i] <- (outmat[2, i] + outmat[2, i - 1])
    }
  }
  cnams <- rep("", length(categories) + 1)
  cnams[1] <- paste("X<=", categories[1], sep = "")
  for (i in 1:(length(categories) - 1)) {
    if (cumulative) 
      cnams[i + 1] <- paste("X<=", categories[i + 1], sep = "")
    else cnams[i + 1] <- paste(categories[i], "<X<=", categories[i + 
                                                                   1], sep = "")
  }
  if (cumulative) 
    cnams[length(categories) + 1] <- paste("all X", sep = "")
  else cnams[length(categories) + 1] <- paste("X>", categories[length(categories)], 
                                              sep = "")
  colnames(outmat) <- cnams
  rownames(outmat) <- c("No", "Prop")
  outmat <- round(outmat, digits = digits)
  outmat
}

