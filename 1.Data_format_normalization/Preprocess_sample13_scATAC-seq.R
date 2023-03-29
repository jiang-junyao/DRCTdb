library(Seurat)
library(Signac)

peaks <- data.table::fread('../Data/scATAC_CSV_file_for_Scanpy/mergedATAC-peaks.csv')

row_granges <- str_replace(peaks[[1]],'-',':')
peaks <- peaks[,-1]
merged_matrix <- Matrix::as.matrix(peaks)
rownames(merged_matrix) <- row_granges

chrom_assay <- CreateChromatinAssay(
  counts = merged_matrix,
  sep = c(":", "-"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)

label <- read.csv("../Data/scATAC_CSV_file_for_Scanpy/mergedATAC-metadata.csv", row.names=1)

na_col <- na.omit(purrr::map_int(1:ncol(label),function(i){
  if (all(is.na(label[[i]]))) {
    return(i)
  }else{
    return(NA)
  }
}))
label <- label[,-na_col]
label_celltype <- label$predicted.id 

sample13 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = label
)

saveRDS(sample13,'sample13_Bone_marrow_ATAC_4k_processed.Rds')
