### Transfer data to seurat object
library(Seurat)
library(tidyverse)
merged_matrix <- ReadMtx(mtx = '../../data/scRNA-seq/sample1/GSE165838_CARE_RNA_counts.mtx.gz',
                         cells = '../../data/scRNA-seq/sample1/GSE165838_CARE_RNA_barcodes.txt.gz',cell.column = 2, skip.cell = 1,
                         features = '../../data/scRNA-seq/sample1/GSE165838_CARE_RNA_features.tsv.gz',feature.column = 2,skip.feature = 1)
metadata <- read.delim("../../data/scRNA-seq/sample1/GSE165838_CARE_RNA_metadata.txt.gz", row.names=1)
obj <- CreateSeuratObject(merged_matrix,meta.data = metadata)

saveRDS(obj,file = '../../data/scRNA-seq/sample1/sample1_heart_scRNA_35k_processed.Rds')


