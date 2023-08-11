library(Signac)
library(Seurat)
library(tidyverse)
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Matrix)
library(glue)
source('../2.Data preprocessing/preprocess_functions.R')
sample1_ATAC <- readRDS('../../data/scATAC-seq/sample1/Rds/sample1_scATAC-seq_80k_processed.Rds')


pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

motif_symbol <- map_chr(1:length(pfm),function(i){
  if (is.null(pfm@listData[[i]]@tags$remap_tf_name)) {
    return(pfm@listData[[i]]@name)
  }else{
    return(pfm@listData[[i]]@tags$remap_tf_name)
  }
})
names(motif_symbol) <- names(pfm)

# add motif information
sample1_ATAC <- AddMotifs(
  object = sample1_ATAC,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

sample1_ATAC <- RunChromVAR(
  object = sample1_ATAC,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
DefaultAssay(sample1_ATAC) <- 'chromvar'


motif_activate_matrix <- sample1_ATAC@assays$chromvar$data
motif_activate_matrix <- as.sparse(motif_activate_matrix)
pseudobulk <- generate_pseudobulk(motif_activate_matrix,group_by = sample1_ATAC$cell_type)


cell_top10_TF_list <- map(pseudobulk, function(x) {
  rownames(pseudobulk)[order(x, decreasing = T)][1:10]
}) 


##match motif
library(motifmatchr)

aCM_gr <- rtracklayer::import.bed('../../data/bed/sample1/aCM.bed.gz')

aCM_motif <- pfm[cell_top10_TF_list$aCM]

motif_ix <- matchMotifs(aCM_motif, aCM_gr, genome = "hg38",out = "positions") 
names(motif_ix) <- map_chr(names(motif_ix),mat2name)

sample <- 'sample1'

for (i in 1:length(motif_ix)) {
  if (!dir.exists(glue('../../data/bed/{sample}'))) {
    dir.create(glue('../../data/bed/{sample}'))
  }
  cell_type <- 'aCM'
  if (!dir.exists(glue('../../data/bed/{sample}/{cell_type}/'))) {
    dir.create(glue('../../data/bed/{sample}/{cell_type}/'))
  }
  filenames <- paste0(glue('../../data/bed/{sample}/{cell_type}/'),names(motif_ix)[i],'.bed.gz')
  rtracklayer::export.bed(object = motif_ix[[i]],con = filenames)
  cat(names(motif_ix)[i],'\n')
}

