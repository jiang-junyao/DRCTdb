library(Seurat)
library(data.table)
library(tidyverse)
##Healthy data
scRNA_Healthy <- readRDS('scRNA-seq/scRNA-Healthy-Hematopoiesis-191120.rds')
meta_data <- as.data.frame(scRNA_Healthy@colData)
meta_data$cell_type <- str_extract(meta_data$BioClassification,'(?<=_).*')
meta_data$cell_type  <- str_extract(meta_data$cell_type,'.*(?<![.\\d])') 

sample16_heathy_RNA <- CreateSeuratObject(
  counts = scRNA_Healthy@assays@.xData$data$counts,
  meta.data = meta_data,
  project = 'Bone_marrow_RNA_Healthy_35k'
)
saveRDS(sample16_heathy_RNA,'scRNA-seq/sample16_Bone_marrow_RNA_Healthy_35k_processed.Rds')
##All data
scRNA_all <- readRDS('scRNA-seq/scRNA-All-Hematopoiesis-MPAL-191120.rds')
meta_data2 <- as.data.frame(scRNA_all@colData)
meta_data2 <- rownames_to_column(meta_data2,var = 'rowname')
meta_data2 <- left_join(meta_data2,meta_data[,c(1,2,9,10)],by = c('Group','Barcode'))
meta_data2 <- as.data.frame(meta_data2)
rownames(meta_data2) <- meta_data2$rowname
get_celltype <- function(i){
  if_else(is.na(meta_data2$BioClassification[i]),
          meta_data2$ProjectClassification[i],
          stringr::str_extract(meta_data2$BioClassification[i],'(?<=_).*')
  )
}
meta_data2$cell_type <- map_chr(1:nrow(meta_data2),get_celltype)

sample16_all_RNA <- CreateSeuratObject(
  counts = scRNA_all@assays@.xData$data$counts,
  meta.data = meta_data2,
  project = 'Bone_marrow_RNA_all_70K'
)


sample16_all_RNA$cell_type <- str_extract(sample16_all_RNA$cell_type,'.*(?<![.\\d])')
saveRDS(sample16_all_RNA,'scRNA-seq/sample16_Bone_marrow_RNA_ALL_56k_processed.Rds')
