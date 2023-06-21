library(Seurat)
library(Signac)
library(tidyverse)
#FCA_GND8046539

load_matrix <- function(x,path = '../../data/scATAC-seq/Sample14/E-MTAB-10570/'){
  mrtrix = paste0(path,x,'_matrix.mtx.gz')
  cells = paste0(path,x,'_barcodes.tsv.gz')
  feature = paste0(path,x,'_features.tsv.gz')
  nrow <- data.table::fread(feature) %>% nrow()

    merged_matrix <- ReadMtx(mtx = mrtrix,
                             features = feature,
                             cells = cells,
                             feature.column = 1)


  colnames(merged_matrix) <- paste0(x,'#',colnames(merged_matrix))
  rna_seurat <- CreateSeuratObject(
    counts = merged_matrix,
    assay = "RNA",
    project = x
  )
  return(rna_seurat)
}

map(all_file,function(x){
  feature = paste0('../../data/scRNA-seq/Sample14/E-MTAB-10551/',x,'_features.tsv.gz')
  data.table::fread(feature) %>% nrow()
})

metdata_list <- map_dfr(metdata_files,data.table::fread)

intersect(unique(metdata_list$sample),unique(str_extract(colnames(sce),'\\w+')))

metdata_list$barcorde <- str_extract(metdata_list$V1,'(?<=\\d_).*')
metdata_list$cells <- paste0(metdata_list$sample,'#', metdata_list$barcorde, '-1') 
cells <-metdata_list$cells %>% intersect(colnames(sce)) 

all_file <- unique(na.omit(str_extract(list.files('../../data/scRNA-seq/Sample14/E-MTAB-10551/'),'.*(?,=_)')))
sce_list <- map(all_file[c(-9)][-21][-60],load_matrix,path = '../../data/scRNA-seq/Sample14/E-MTAB-10551/',.progress = T)
sce <- purrr::reduce(sce_list,merge)


sample14_rna <- sce[,cells]
sample14_metadada <- metdata_list[metdata_list]

