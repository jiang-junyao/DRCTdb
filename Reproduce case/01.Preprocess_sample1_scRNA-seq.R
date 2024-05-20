### Transfer data to seurat object
library(Seurat)
library(tidyverse)
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165838
merged_matrix <- ReadMtx(mtx = 'GSE165838_CARE_RNA_counts.mtx.gz',
                         cells = 'GSE165838_CARE_RNA_barcodes.txt.gz',cell.column = 2, skip.cell = 1,
                         features = 'GSE165838_CARE_RNA_features.tsv.gz',feature.column = 2,skip.feature = 1)
metadata <- read.delim("GSE165838_CARE_RNA_metadata.txt.gz", row.names=1)
obj <- CreateSeuratObject(merged_matrix,meta.data = metadata)

saveRDS(obj,file = '../../data/scRNA-seq/sample1/sample1_heart_scRNA_35k_processed.Rds')


