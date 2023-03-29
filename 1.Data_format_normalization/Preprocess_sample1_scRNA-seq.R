### Transfer data to seurat object
library(Seurat)
library(tidyverse)
merged_matrix <- ReadMtx(mtx = './rawdata/GSE165838_CARE_RNA_counts.mtx.gz',
                         cells = './rawdata/GSE165838_CARE_RNA_barcodes.txt.gz',cell.column = 2, skip.cell = 1,
                         features = './rawdata/GSE165838_CARE_RNA_features.tsv.gz',feature.column = 2,skip.feature = 1)
metadata <- read.delim("./rawdata/GSE165838_CARE_RNA_metadata.txt.gz", row.names=1)
obj <- CreateSeuratObject(merged_matrix,meta.data = metadata)


obj <- SCTransform(obj) %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:30,umap.method = 'umap-learn',metric = 'correlation',verbose = FALSE) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE) %>%
  FindClusters(verbose = FALSE)
saveRDS(obj,file = 'Rds/heart_processed.Rds')


