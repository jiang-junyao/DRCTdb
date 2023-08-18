library(tidyverse)
library(Seurat)
sce_list <- map(list.files('./scRNA-seq/',pattern = '*.rds',full.names = T),function(x){
    sce <- readRDS(x)
    colnames(sce) <- paste0(colnames(sce),'#',sce$orig.ident)
    DefaultAssay(sce) <- 'RNA'
    sce <- DietSeurat(sce,assays = 'RNA')
    cat(x,'finished\n')
    return(sce)
})

sce <- reduce(sce_list,merge)
sce$cell_type <- sce$CellType
saveRDS(sce,'sample27_processed_280kRNA.Rds')
