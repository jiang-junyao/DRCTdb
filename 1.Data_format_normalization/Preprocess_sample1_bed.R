library(Signac)
library(Seurat)
library(tidyverse)
library(Matrix)


sample1_ATAC <- readRDS('../../data/scATAC-seq/sample1/Rds/sample1_scATAC-seq_80k_processed.Rds')

df <- sample1_ATAC@assays$peaks@counts
df <- df[str_starts(rownames(df),'chr'),]

##subset columns, if your peak is less than 1M, your does not need to use multithreads
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
pseudobulk <- generate_pseudobulk(df,group_by = sample1_ATAC$cell_type)

# pseudobulk <- generate_pseudobulk(df,group_by = na.omit(sample1_ATAC$cell_type))
# identical(rownames(pseudobulk),rownames(df))
# system.time(generate_pseudobulk(df,group_by = na.omit(sample1_ATAC$cell_type)))
# 
# pseudobulk2 <- generate_pseudobulk(df,group_by = na.omit(sample1_ATAC$cell_type),ncore = 10)
# identical(rownames(pseudobulk),rownames(df))
# identical(rownames(pseudobulk),rownames(pseudobulk2))
# system.time(generate_pseudobulk(df,group_by = na.omit(sample1_ATAC$cell_type),ncore = 10))

cell_gr <- separate(as.data.frame(rownames(pseudobulk)),col = everything(),sep = '-',into = c('seqnames','start','end'))

cell_gr_list <- purrr::map(pseudobulk,function(x){
  threshold <- str_extract(colnames(as.data.frame(catable(x)))[5],'\\d+') %>% as.numeric()
  gr <- cell_gr[which(as.vector(x) >= threshold),] %>% GenomicRanges::makeGRangesFromDataFrame()
  return(gr)
})

for (i in 1:length(cell_gr_list)) {
  filenames <- paste0('../../data/bed/sample1/',names(cell_gr_list)[i],'.bed.gz')
  rtracklayer::export.bed(object = cell_gr_list[[i]],con = filenames)
  cat(names(cell_gr_list)[i],'\n')
}

