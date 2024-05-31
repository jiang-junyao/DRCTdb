### Transfer data to seurat object
library(Seurat)
library(tidyverse)
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165838
merged_matrix <- ReadMtx(mtx = 'GSE165838_CARE_RNA_counts.mtx.gz',
                         cells = 'GSE165838_CARE_RNA_barcodes.txt.gz',cell.column = 2, skip.cell = 1,
                         features = 'GSE165838_CARE_RNA_features.tsv.gz',feature.column = 2,skip.feature = 1)
metadata <- read.delim("GSE165838_CARE_RNA_metadata.txt.gz", row.names=1)
obj <- CreateSeuratObject(merged_matrix,meta.data = metadata)

saveRDS(obj,file = 'sample1_heart_scRNA_35k_processed.Rds')


library(Seurat)
library(Signac)
merged_matrix <- ReadMtx('GSE165837_CARE_ATAC_merged_matrix.mtx.gz',
                         cells = 'GSE165837_CARE_ATAC_merged_barcodes.txt.gz',
                         features = 'GSE165837_CARE_ATAC_merged_features.txt.gz',
                         feature.column = 1)

merged_matrix <- merged_matrix[str_starts(rownames(merged_matrix),'chr'),]


label <- read.delim("GSE165837_CARE_ATAC_merged_umap_cluster_labels.tsv.gz", row.names=1)

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
saveRDS(sample1,file = 'sample1_scATAC-seq_80k_processed.Rds')
dior::seurat_write_h5(sample1, file='sample1_scATAC-seq_80k_processed.h5',assay.name = DefaultAssay(sample1),save.graphs = F)
