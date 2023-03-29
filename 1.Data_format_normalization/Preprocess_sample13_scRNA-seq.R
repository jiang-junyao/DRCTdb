library(Seurat)
library(tidyverse)
merged_matrix <- ReadMtx(mtx = 'MergedAllSamples_annotated/MergedAllSamples_annotated.mtx.gz',
                         cells = './MergedAllSamples_annotated/MergedAllSamples_annotated_metadata.txt.gz',cell.column = 1, skip.cell = 1,
                         features = 'MergedAllSamples_annotated/MergedAllSamples_annotated_features.txt.gz',feature.column = 1,skip.feature = 1)
metadata <- read.delim("./MergedAllSamples_annotated/MergedAllSamples_annotated_metadata.txt.gz", row.names=1)


obj <- CreateSeuratObject(merged_matrix,meta.data = metadata,project = 'Bone marrow 8k')
obj$cell_type <- str_extract(obj$manual_annotation,'.*(?<!\\s\\d)') 
saveRDS(obj,'sample13_Bone_marrow_RNA_4k_processed.Rds')
##
# sce <- obj %>% SCTransform() %>% RunPCA(verbose = FALSE) 
# sce <- sce %>% RunUMAP(dims = 1:30,umap.method = 'umap-learn',metric = 'correlation',verbose = FALSE) 
# sce$cell_type <- str_extract(sce$manual_annotation,'.*(?<!\\s\\d)') 
# DimPlot(sce,group.by = 'cell_type')
