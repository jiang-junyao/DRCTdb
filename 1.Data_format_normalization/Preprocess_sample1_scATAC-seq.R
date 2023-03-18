### Transfer data to seurat object
library(Seurat)

merged_matrix <- ReadMtx('../ignore/GSE165837_CARE_ATAC_merged_matrix.mtx.gz',
                         cells = '../ignore/GSE165837_CARE_ATAC_merged_barcodes.txt.gz',
                         features = '../ignore/GSE165837_CARE_ATAC_merged_features.txt.gz',
                         feature.column = 1)
label <- read.delim("E:/DRCTdb/ignore/GSE165837_CARE_ATAC_merged_umap_cluster_labels.tsv.gz", row.names=1)
obj <- CreateSeuratObject(merged_matrix,meta.data = label)
saveRDS(merged_matrix,file = '../data/scATAC-seq/sample1/merged_matrix.Rds')
