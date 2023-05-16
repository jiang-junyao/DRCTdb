library(ArchR)
addArchRThreads(threads = 1) 
addArchRGenome("hg38")
sample4 <- loadArchRProject('../../data/scATAC-seq/Sample4/',showLogo = F)

peak_matrix <- getMatrixFromProject(
  ArchRProj = sample4,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE
)
saveRDS(peak_matrix,file = '../../data/scATAC-seq/Sample4/Sample4_peak_matrix.Rds')
library(tidyverse)
library(Matrix)
library(glue)
source('preprocess_functions.R')
sample4_ATAC <- readRDS('../../data/scATAC-seq/Sample4/Sample4_peak_matrix.Rds')
sparse_mtx <- sample4_ATAC@assays@data$PeakMatrix
rownames(sparse_mtx) <- paste(as.data.frame(sample4_ATAC@rowRanges)[[1]],as.data.frame(sample4_ATAC@rowRanges)[[2]],as.data.frame(sample4_ATAC@rowRanges)[[3]],sep = '-')
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sparse_mtx),subset_peaks)),]

pseudobulk <- generate_pseudobulk(sparse_mtx,group_by = sample4_ATAC$Celltype)

cell_gr <- tidyr::separate(as.data.frame(rownames(pseudobulk)),col = everything(),sep = '-',into = c('seqnames','start','end'))
cell_gr_list <- map(colnames(pseudobulk),get_cell_gr,seurat_Obj = sample4_ATAC)
names(cell_gr_list) <- colnames(pseudobulk)

sample <- 'sample4'

for (i in 1:length(cell_gr_list)) {
  if (!dir.exists(glue('../../data/bed/{sample}'))) {
    dir.create(glue('../../data/bed/{sample}'))
  }
  filenames <- paste0(glue('../../data/bed/{sample}/'),names(cell_gr_list)[i],'.bed.gz')
  rtracklayer::export.bed(object = cell_gr_list[[i]],con = filenames)
  cat(names(cell_gr_list)[i],'\n')
}

