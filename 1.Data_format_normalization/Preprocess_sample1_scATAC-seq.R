### Transfer data to seurat object
library(Seurat)
library(Signac)
merged_matrix <- ReadMtx('../../data/scATAC-seq/sample1/GSE165837_CARE_ATAC_merged_matrix.mtx.gz',
                         cells = '../../data/scATAC-seq/sample1/GSE165837_CARE_ATAC_merged_barcodes.txt.gz',
                         features = '../../data/scATAC-seq/sample1/GSE165837_CARE_ATAC_merged_features.txt.gz',
                         feature.column = 1)
label <- read.delim("../../data/scATAC-seq/sample1/GSE165837_CARE_ATAC_merged_umap_cluster_labels.tsv.gz", row.names=1)

chrom_assay <- CreateChromatinAssay(
  counts = merged_matrix,
  sep = c(":", "-"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)
sample1 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = label
)
sample1$cell_type <- sample1$Cluster
saveRDS(sample1,file = '../../data/scATAC-seq/sample1/Rds/sample1_seurat_obj.Rds')

