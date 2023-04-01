library(Seurat)
library(Signac)
library(data.table)
library(tidyverse)
##Healthy data
scATAC_Healthy <- readRDS('scATAC-seq/Github/scATAC-Healthy-Hematopoiesis-191120.rds')
meta_data <- as.data.frame(scATAC_Healthy@colData)
meta_data$cell_type <- stringr::str_extract(meta_data$BioClassification,'(?<=_).*')

rtracklayer::export.bed(scATAC_Healthy@rowRanges,'scATAC-seq/Github/scATAC_Healthy_hg19_cord.bed')

hg38 <- fread('scATAC-seq/Github/scATAC_Healthy_hg38_lift_cord.bed') |>  unique(by="V4")



sparce_mtx <- scATAC_Healthy@assays@.xData$data$counts[-which(!(scATAC_Healthy@rowRanges$name %in% hg38$V4)),]
#Remove missing 
dim(sparce_mtx)
colnames(sparce_mtx) <- rownames(meta_data)
rownames(sparce_mtx) <- paste0(hg38$V1,':',hg38$V2,'-',hg38$V3)


chrom_assay <- CreateChromatinAssay(
  counts = sparce_mtx,
  sep = c(":", "-"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)
sample16 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = meta_data
)
saveRDS(sample16,'scATAC-seq/sample16_Bone_marrow_ATAC_Healthy_35k_processed.Rds')



####Healthy data-----
scATAC_all <- readRDS('scATAC-seq/Github/scATAC-All-Hematopoiesis-MPAL-191120.rds')
scATAC_all
meta_data2 <- as.data.frame(scATAC_all@colData)
meta_data2 <- left_join(meta_data2,meta_data[,c(3,4,9)],by = c('Group','Barcode'))

get_celltype <- function(i){
  if_else(is.na(meta_data2$BioClassification[i]),
          meta_data2$ProjectClassification[i],
          stringr::str_extract(meta_data2$BioClassification[i],'(?<=_).*')
  )
}
meta_data2$cell_type <- map_chr(1:nrow(meta_data2),get_celltype)


rtracklayer::export.bed(scATAC_all@rowRanges,'scATAC-seq/Github/scATAC_all_hg19_cord.bed')

hg38_all <- fread('scATAC-seq/Github/scATAC_all_hg38_lift_cord.bed') |>  unique(by="V4")
sparce_mtx2 <- scATAC_all@assays@.xData$data$counts[-which(!(scATAC_all@rowRanges$name %in% hg38_all$V4)),]
#Remove missing 
dim(sparce_mtx2)
colnames(sparce_mtx2) <- rownames(meta_data2)
rownames(sparce_mtx2) <- paste0(hg38_all$V1,':',hg38_all$V2,'-',hg38_all$V3)

chrom_assay2 <- CreateChromatinAssay(
  counts = sparce_mtx2,
  sep = c(":", "-"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)
sample16_2 <- CreateSeuratObject(
  counts = chrom_assay2,
  assay = "peaks",
  meta.data = meta_data2
)
saveRDS(sample16_2,'scATAC-seq/sample16_Bone_marrow_ATAC_all_70k_processed.Rds')
