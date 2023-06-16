library(Seurat)
library(Signac)
library(tidyverse)
#FCA_GND8046539

load_matrix <- function(x){
  mrtrix = paste0('../../../../Desktop/Sample14/E-MTAB-10570/',x,'_matrix.mtx.gz')
  feature = paste0('../../../../Desktop/Sample14/E-MTAB-10570/',x,'_barcodes.tsv.gz')
  cells = paste0('../../../../Desktop/Sample14/E-MTAB-10570/',x,'_features.tsv.gz')
  merged_matrix <- ReadMtx(mtx = mrtrix,
                           features = feature,
                           cells = cells,
                           feature.column = 1)
  merged_matrix <- Matrix::t(merged_matrix)
  chrom_assay <- CreateChromatinAssay(
    counts = merged_matrix,
    sep = c(":", "-"),
    genome = 'hg38',
    min.cells = 10,
    min.features = 200
  )
  return(chrom_assay)
}

all_file <- unique(na.omit(str_extract(list.files('../../../../Desktop/Sample14/E-MTAB-10570/'),'.*(?=_)')))

chrom_assay_list <- map(all_file,load_matrix)

chrom_assay <- reduce(chrom_assay_list,merge)

sample14 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

saveRDS(sample13,'../../../../Desktop/Sample14/peak_matrix.Rds')
