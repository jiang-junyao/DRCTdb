### Transfer data to seurat object
library(Seurat)
library(tidyverse)
sce_list <- map(list.files('../../data/scRNA-seq/Sample28/drive-download-20231030T100413Z-001/',full.names = T)[1:5],function(x){
    sce <- readRDS(x)
    sce$cell_type <- str_extract(basename(x),'\\w+')
    return(sce)
},.progress = T) 

sample28_RNA <- reduce(sce_list,merge)

saveRDS(sample28_RNA,'../../data/scRNA-seq/Sample28/sample28_scRNA_46K_processed.Rds')
library(reticulate)
use_condaenv('scanpy')
loompy <- reticulate::import('loompy')
sceasy::convertFormat(sample28_RNA, from="seurat", to="anndata",assay = 'RNA',
                      outFile='../../data/scRNA-seq/Sample28/sample28_scRNA_46K_processed.h5ad')
