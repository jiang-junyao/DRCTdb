library(Signac)
library(Seurat)
library(tidyverse)
library(Matrix)
library(glue)
source('preprocess_functions.R')
sample11_ATAC <- readRDS('../../data/scATAC-seq/Sample11/Sample11_peak_matrix.Rds')
sample11_ATAC@colData$barcode <- rownames(sample11_ATAC@colData)
sample11_ATAC@colData$barcode <- str_extract(sample11_ATAC@colData$barcode,'(?<=#).*')
  
coldata <- as.data.frame(sample11_ATAC@colData)
coldata$raw_barcode <- rownames(coldata)
sparse_mtx <- sample11_ATAC@assays@data$PeakMatrix
rownames(sparse_mtx) <- paste(as.data.frame(sample11_ATAC@rowRanges)[[1]],as.data.frame(sample11_ATAC@rowRanges)[[2]],as.data.frame(sample11_ATAC@rowRanges)[[3]],sep = '-')
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sparse_mtx),subset_peaks)),]

metadata <- readRDS('../../data/scATAC-seq/Sample11/meta.rds')
metadata$barcode <- rownames(metadata)
metadata$barcode <- str_replace(metadata$barcode,'h_Donor','hr_Donor')

coldata <- coldata %>% left_join(metadata,by = 'barcode')
coldata$cell_type <- coldata$pairedLabel

coldata2 <- coldata[which(!is.na(coldata$cell_type)),]
sparse_mtx <- sparse_mtx[,which(!is.na(coldata$cell_type))]
identical(colnames(sparse_mtx),coldata2$raw_barcode)

pseudobulk <- generate_pseudobulk(sparse_mtx,group_by = coldata2$cell_type)
cell_gr <- tidyr::separate(as.data.frame(rownames(pseudobulk)),col = everything(),sep = '-',into = c('seqnames','start','end'))
cell_gr_list <- map(colnames(pseudobulk),get_cell_gr,seurat_Obj = coldata2)
names(cell_gr_list) <- colnames(pseudobulk)

sample <- 'sample11'

for (i in 1:length(cell_gr_list)) {
  if (!dir.exists(glue('../../data/bed/{sample}'))) {
    dir.create(glue('../../data/bed/{sample}'))
  }
  filenames <- paste0(glue('../../data/bed/{sample}/'),names(cell_gr_list)[i],'.bed.gz')
  rtracklayer::export.bed(object = cell_gr_list[[i]],con = filenames)
  cat(names(cell_gr_list)[i],'\n')
}

