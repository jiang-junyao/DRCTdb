### Transfer data to seurat object
library(Seurat)
library(Signac)
label <- data.table::fread("../../data/scATAC-seq/Sample23/multiome/GSE200044_multiome_beta_annotation.csv.gz") |> as.data.frame()
rownames(label) <- label$X
label <- label[,c(-1,-2)]
#scRNA-seq
merged_matrix_RNA <- ReadMtx('../../data/scATAC-seq/Sample23/multiome/beta_RNA.mtx',
                         features = '../../data/scATAC-seq/Sample23/multiome/beta_RNA.barcodes',
                         cells = '../../data/scATAC-seq/Sample23/multiome/beta_RNA.genes',
                         feature.column = 1)
merged_matrix_RNA <- Matrix::t(merged_matrix_RNA)
obj <- CreateSeuratObject(merged_matrix_RNA,meta.data = label)

saveRDS(merged_matrix_RNA,'../../data/scRNA-seq/Sample23/sample23_scRNAseq.Rds')

#scATAC-seq
merged_matrix_ATAC <- ReadMtx('../../data/scATAC-seq/Sample23/multiome/beta_ATAC_500bp.mtx',
                         features = '../../data/scATAC-seq/Sample23/multiome/beta_ATAC_500bp.barcodes',
                         cells = '../../data/scATAC-seq/Sample23/multiome/beta_ATAC_500bp.regions',
                         feature.column = 1)
merged_matrix_ATAC <- Matrix::t(merged_matrix_ATAC)


chrom_assay <- CreateChromatinAssay(
  counts = merged_matrix_ATAC,
  sep = c(":", "-"),
  genome = 'hg19',
  min.cells = 10,
  min.features = 200
)
sample23 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = label
)
sample23$cell_type <- sample23$subtype
saveRDS(sample23,file = '../../data/scATAC-seq/sample23/sample23_scATAC-seq.Rds')
