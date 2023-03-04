library(Seurat)
merged_matrix <- ReadMtx('../data/scATAC-seq/sample1/matrix.mtx.gz',
                         cells = '../data/scATAC-seq/sample1/barcodes.tsv.gz',
                         features = '../data/scATAC-seq/sample1/features.tsv.gz',
                         feature.column = 1)
save(merged_matrix,file = '../data/scATAC-seq/sample1/merged_matrix.Rds')
