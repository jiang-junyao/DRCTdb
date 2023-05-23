library(Seurat)
df <- readRDS('human_6_8_12and19_merged_final_cleaned.rds')
mat <- df@assays$RNA@counts
meta <- df@meta.data
sample4 <- CreateSeuratObject(mat,assay = 'RNA',project = 'sample4')
sample4 <- AddMetaData(sample4,metadata = meta)
sample4 <- SCTransform(sample4) |> RunPCA(verbose = FALSE) |>
  RunUMAP(dims = 1:30,umap.method = 'umap-learn',metric = 'correlation',verbose = FALSE)
sample4$cell_type <- sample4$new_manual_annotation
DimPlot(sample4,group.by = 'new_manual_annotation')
saveRDS(sample4,file = 'sample4_scRNA-seq.Rds')