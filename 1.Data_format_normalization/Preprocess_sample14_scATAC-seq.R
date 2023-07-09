library(Seurat)
library(Signac)
library(tidyverse)
#FCA_GND8046539

load_matrix <- function(x,path = '../../data/scATAC-seq/Sample14/E-MTAB-10570/'){
  mrtrix = paste0(path,x,'_matrix.mtx.gz')
  feature = paste0(path,x,'_barcodes.tsv.gz')
  cells = paste0(path,x,'_features.tsv.gz')
  merged_matrix <- ReadMtx(mtx = mrtrix,
                           features = feature,
                           cells = cells,
                           feature.column = 1)
  merged_matrix <- Matrix::t(merged_matrix)
  colnames(merged_matrix) <- paste0(x,'#',colnames(merged_matrix))
  chrom_assay <- CreateChromatinAssay(
    counts = merged_matrix,
    sep = c(":", "-"),
    genome = 'hg38',
    min.cells = 10,
    min.features = 200
  )
  return(chrom_assay)
}

all_file <- unique(na.omit(str_extract(list.files('../../data/scATAC-seq/Sample14/E-MTAB-10570/'),'.*(?=_)')))

chrom_assay_list <- map(all_file,load_matrix)

chrom_assay <- purrr::reduce(chrom_assay_list,merge)

sample14 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)
metadata <- data.table::fread('../../data/scATAC-seq/Sample14/atac_metadata.csv',header = T)

metadata2 <- data.frame(
  barcode = paste0(metadata$sample,'#',str_extract(metadata$V1,'[A-Z]+'),'-1'),
  cell_type = metadata$cell_type
)
rownames(metadata2) <- metadata2$barcode


cells <- intersect(colnames(sample14),metadata2$barcode)
sample14_2 <- subset(sample14,cells = cells)
metadata2 <- metadata2[cells,]

identical(colnames(sample14_2),rownames(metadata2))

sample14_2 <- AddMetaData(sample14_2,metadata = metadata2)
saveRDS(sample14_2,'../../data/scATAC-seq/Sample14/sample14_processed_95k_scATAC.Rds')
##multiomics---
all_file <- unique(na.omit(str_extract(list.files('../../data/scATAC-seq/Sample14/multiome/E-MTAB-11708/'),'.*(?=_)')))
load_matrix2 <- function(x){
  mrtrix = paste0('../../data/scATAC-seq/Sample14/multiome/E-MTAB-11708/',x,'_matrix.mtx.gz')
  cells= paste0('../../data/scATAC-seq/Sample14/multiome/E-MTAB-11708/',x,'_barcodes.tsv.gz')
  feature = paste0('../../data/scATAC-seq/Sample14/multiome/E-MTAB-11708/',x,'_features.tsv.gz')
  merged_matrix <- ReadMtx(mtx = mrtrix,
                           features = feature,
                           cells = cells,
                           feature.column = 1)
  colnames(merged_matrix) <- paste0(x,'#',colnames(merged_matrix))
  
  
  RNA <- merged_matrix[str_starts(rownames(merged_matrix),'ENSG'),]
  ATAC <- merged_matrix[!str_starts(rownames(merged_matrix),'ENSG'),]
  chrom_assay <- CreateChromatinAssay(
    counts = ATAC,
    sep = c(":", "-"),
    genome = 'hg38',
    min.cells = 0,
    min.features = 0
  )
  rna_seurat <- CreateSeuratObject(
    counts = RNA,
    assay = "RNA"
  )
  rna_seurat[['"ATAC"']] <- chrom_assay
  return(rna_seurat)
}

sample14_multiome_list <- map(all_file,load_matrix2)
sample14_multiome <- purrr::reduce(sample14_multiome_list,merge)

metdata_files <- list.files('../../data/scATAC-seq/Sample14/metadata/',full.names = T)
metdata_list <- map_dfr(metdata_files,data.table::fread)

intersect(unique(metdata_list$sample),unique(str_extract(colnames(sample14_multiome),'\\w+')))
saveRDS(sample14_multiome,'../../data/scATAC-seq/Sample14/sample14_processed_50k_multiome.Rds')





