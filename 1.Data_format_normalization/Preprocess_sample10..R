library(Seurat)
library(Signac)
merged_matrix <- ReadMtx('sample10_fetal_lung/sample10_fetal_lung.mtx',
                         cells = 'sample10_fetal_lung/sample10_fetal_lung_features.txt.gz',
                         features = 'sample10_fetal_lung/sample10_fetal_lung_metadata.txt.gz',
                         feature.column = 1)
