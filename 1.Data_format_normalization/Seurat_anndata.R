library(Seurat)
library(Signac)
library(sceasy)
library(tidyverse)
library(reticulate)
use_condaenv('scanpy')
loompy <- reticulate::import('loompy')
source('../2.Data preprocessing/preprocess_functions.R')
#sample1----
sample1_ATAC <- readRDS('../../data/scATAC-seq/sample1/Rds/sample1_scATAC-seq_80k_processed.Rds')
sceasy::convertFormat(sample1_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/sample1/Rds/sample1_scATAC-seq_80k_processed.h5ad')
#sample4----
sample3_mtx <- readRDS('../../data/scATAC-seq/Sample3/sample3_peak_matrix.Rds')
sparse_mtx <- sample3_mtx@assays@data$PeakMatrix
rownames(sparse_mtx) <- paste(as.data.frame(sample3_mtx@rowRanges)[[1]],as.data.frame(sample3_mtx@rowRanges)[[2]],as.data.frame(sample3_mtx@rowRanges)[[3]],sep = '-')
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sparse_mtx),subset_peaks)),]

chrom_assay <- CreateChromatinAssay(
    counts = sparse_mtx,
    sep = c("-", "-"),
    genome = 'hg38',
    min.cells = 10,
    min.features = 200
)
sample3_ATAC <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    metadata = as.data.frame(sample3_mtx@colData)
)
identical(colnames(sparse_mtx),rownames(sample3_mtx@colData))
sample3_ATAC <- AddMetaData(sample3_ATAC,metadata = as.data.frame(sample3_mtx@colData))
sample3_ATAC$cell_type <- sample3_ATAC$Sample
saveRDS(sample3_ATAC,file = '../../data/scATAC-seq/Sample3/sample3_scATAC-seq_756k_processed.Rds')
sceasy::convertFormat(sample3_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/Sample3/sample3_scATAC-seq_756k_processed.h5ad')

#sample4----
sample4_mtx <- readRDS('../../data/scATAC-seq/Sample4/Sample4_peak_matrix.Rds')
sparse_mtx <- sample4_mtx@assays@data$PeakMatrix
rownames(sparse_mtx) <- paste(as.data.frame(sample4_mtx@rowRanges)[[1]],as.data.frame(sample4_mtx@rowRanges)[[2]],as.data.frame(sample4_mtx@rowRanges)[[3]],sep = '-')
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sparse_mtx),subset_peaks)),]
sample4_mtx$cell_type <- sample4_mtx$Celltype
chrom_assay <- CreateChromatinAssay(
  counts = sparse_mtx,
  sep = c("-", "-"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)
sample4_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)
identical(colnames(sparse_mtx),rownames(sample4_mtx@colData))
sample4_ATAC <- AddMetaData(sample4_ATAC,metadata = as.data.frame(sample4_mtx@colData))
saveRDS(sample4_ATAC,file = '../../data/scATAC-seq/Sample4/sample4_scATAC-seq_30k_processed.Rds')
sceasy::convertFormat(sample4_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/sample4/sample4_scATAC-seq_30k_processed.h5ad')
#sample5----
Sample5_rds <- list.files('../../data/scATAC-seq/Sample5/',pattern = 'RDS',full.names = T)
sample5_sce <- map(Sample5_rds,function(x){
    sce <- UpdateSeuratObject(readRDS(x))
    return(sce)
    })

for (i in seq_along(sample5_sce)) {
    filename = str_extract(Sample5_rds,'.*(?=.RDS)')[i]
    sceasy::convertFormat(sample5_sce[[i]], from="seurat", to="anndata",assay = 'peaks', outFile=paste0(filename,'.h5ad')) 
}




#sample7----
sample7_mtx<- readRDS('../../data/scATAC-seq/Sample7/Sample7_peak_matrix.Rds')
sparse_mtx <- sample7_mtx@assays@data$PeakMatrix
sample7_mtx$cell_type <- sample7_mtx$Sample
rownames(sparse_mtx) <- paste(as.data.frame(sample7_mtx@rowRanges)[[1]],as.data.frame(sample7_mtx@rowRanges)[[2]],as.data.frame(sample7_mtx@rowRanges)[[3]],sep = '-')
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sparse_mtx),subset_peaks)),]
chrom_assay <- CreateChromatinAssay(
  counts = sparse_mtx,
  sep = c("-", "-"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)
sample7_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = as.data.frame(sample7_mtx@colData)
)
saveRDS(sample7_ATAC,file = '../../data/scATAC-seq/Sample7/sample7_scATAC-seq_50k_processed.Rds')
sceasy::convertFormat(sample7_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/sample7/sample7_scATAC-seq_50k_processed.h5ad')
#sample8----
sample8_ATAC <- readRDS('../../data/scATAC-seq/Sample8/sample8_peak_matrix.Rds')
sparse_mtx <- sample8_ATAC@assays@data$PeakMatrix
rownames(sparse_mtx) <- paste(as.data.frame(sample8_ATAC@rowRanges)[[1]],as.data.frame(sample8_ATAC@rowRanges)[[2]],as.data.frame(sample8_ATAC@rowRanges)[[3]],sep = '-')

chrom_assay <- CreateChromatinAssay(
  counts = sparse_mtx,
  sep = c("-", "-"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)
sample8_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = as.data.frame(sample8_ATAC@colData)
)
saveRDS(sample8_ATAC,file = '../../data/scATAC-seq/Sample8/sample8_scATAC-seq_100k_processed.Rds')
sceasy::convertFormat(sample8_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/sample8/sample8_scATAC-seq_100k_processed.h5ad')
#sample10----
sample10_ATAC <- readRDS('../../data/scATAC-seq/Sample10/sample10_peak_matrix.Rds')
sparse_mtx <- sample10_ATAC@assays@data$PeakMatrix
rownames(sparse_mtx) <- paste(as.data.frame(sample10_ATAC@rowRanges)[[1]],as.data.frame(sample10_ATAC@rowRanges)[[2]],as.data.frame(sample10_ATAC@rowRanges)[[3]],sep = '-')
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sparse_mtx),subset_peaks)),]

chrom_assay <- CreateChromatinAssay(
  counts = sparse_mtx,
  sep = c("-", "-"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)
sample10_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = as.data.frame(sample10_ATAC@colData)
)
sample10_ATAC$cell_type <- str_replace(sample10_ATAC$cell_type,'\\+','_high')
sample10_ATAC$cell_type <- str_replace_all(sample10_ATAC$cell_type,'/','_')
sample10_ATAC$cell_type <- str_replace_all(sample10_ATAC$cell_type,' ','_')
sample10_ATAC$cell_type <- str_replace_all(sample10_ATAC$cell_type,'[.]','')
saveRDS(sample10_ATAC,file = '../../data/scATAC-seq/Sample10/sample10_scATAC-seq_95k_processed.Rds')
sceasy::convertFormat(sample10_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/sample10/sample10_scATAC-seq_95k_processed.h5ad')
#sample11----
sample11_mtx <- readRDS('../../data/scATAC-seq/Sample11/Sample11_peak_matrix.Rds')

sparse_mtx <- sample11_mtx@assays@data$PeakMatrix
rownames(sparse_mtx) <- paste(as.data.frame(sample11_mtx@rowRanges)[[1]],as.data.frame(sample11_mtx@rowRanges)[[2]],as.data.frame(sample11_mtx@rowRanges)[[3]],sep = '-')
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sparse_mtx),subset_peaks)),]


chrom_assay <- CreateChromatinAssay(
  counts = sparse_mtx,
  sep = c("-", "-"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)
sample11_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = as.data.frame(sample11_mtx@colData)
)
sample11_ATAC <- sample11_ATAC[,which(!is.na(sample11_ATAC$cell_type))]

saveRDS(sample11_ATAC,file = '../../data/scATAC-seq/Sample11/sample11_scATAC-seq_17k_processed.Rds')
sceasy::convertFormat(sample11_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/sample11/sample11_scATAC-seq_17k_processed.h5ad')

#sample12----
sample12_ATAC <- readRDS('../../data/scATAC-seq/Sample12/Sample12_ATAC-seq.Rds')
sample12_ATAC$cell_type <- str_replace_all(sample12_ATAC$cell_type,' ','_')
sample12_ATAC$cell_type <- str_replace_all(sample12_ATAC$cell_type,'/','')
sceasy::convertFormat(sample12_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/sample12/sample12_scATAC-seq_104k_processed.h5ad')
#sample13----
sample13_RNA <- readRDS('../../data/scRNA-seq/Sample13/sample13_Bone_marrow_8k.Rds')
sample13_ATAC <- readRDS('../../data/scATAC-seq/Sample13/sample13_Bone_marrow_ATAC_4k_processed.Rds.rds')
sample13_ATAC$cell_type <- sample13_ATAC$predicted.id 
saveRDS(sample13_ATAC,'../../data/scATAC-seq/Sample13/sample13_Bone_marrow_ATAC_4k_processed.Rds')
sceasy::convertFormat(sample13_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/sample13/sample13_scATAC-seq_4k_processed.h5ad')


#sample14----
sample14_ATAC <- readRDS('../../data/scATAC-seq/sample14/sample14_processed_95k_scATAC.Rds')
sceasy::convertFormat(sample14_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/sample14/sample14_scATAC-seq_95k_processed.h5ad')

#sample16----
sample16_ATAC <- readRDS('../../data/scATAC-seq/sample16/sample16_Bone_marrow_ATAC_Healthy_35k_processed.Rds')
sparse_mtx <- sample16_ATAC@assays$peaks@counts
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sparse_mtx),subset_peaks)),]
chrom_assay <- CreateChromatinAssay(
  counts = sparse_mtx,
  sep = c("-", "-"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)
sample16_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = as.data.frame(sample16_ATAC@meta.data)
)

sceasy::convertFormat(sample16_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/sample16/sample16_scATAC-seq_35k_processed.h5ad')
#sample23----
sample23_ATAC <- readRDS('../../data/scATAC-seq/sample23/islet_multiome_seurat_singac.Rds')
DefaultAssay(sample23_ATAC) <- 'ATAC'
sparse_mtx <- sample23_ATAC@assays$ATAC@counts
chrom_assay <- CreateChromatinAssay(
  counts = sparse_mtx,
  sep = c("-", "-"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)
sample23_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = as.data.frame(sample23_ATAC@meta.data)
)

sceasy::convertFormat(sample23_ATAC, from="seurat", to="anndata",assay = DefaultAssay(sample23_ATAC),
                      outFile='../../data/scATAC-seq/sample23/sample23_scATAC-seq_95k_processed.h5ad')
#sample24----
sample24_ATAC <- readRDS('../../data/scATAC-seq/Sample24/sample24_scATAC-seq_100k_healthy_processed.Rds')
sceasy::convertFormat(sample24_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/sample24/sample24_scATAC-seq_52k_processed.h5ad')
#sample24----
sample26_mtx <- readRDS('../../data/scATAC-seq/Sample26/Sample26_peak_matrix.Rds')
sparse_mtx <- sample26_mtx@assays@data$PeakMatrix
rownames(sparse_mtx) <- paste(as.data.frame(sample26_mtx@rowRanges)[[1]],as.data.frame(sample26_mtx@rowRanges)[[2]],as.data.frame(sample26_mtx@rowRanges)[[3]],sep = '-')
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sparse_mtx),subset_peaks)),]
chrom_assay <- CreateChromatinAssay(
    counts = sparse_mtx,
    sep = c("-", "-"),
    genome = 'hg19',
    min.cells = 10,
    min.features = 200
)
sample26_ATAC <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = as.data.frame(sample26_mtx@colData)
)
saveRDS(sample26_ATAC,file = '../../data/scATAC-seq/Sample26')
sceasy::convertFormat(sample26_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/Sample26/Sample26_scATAC-seq_229k_processed.h5ad')

sceasy::convertFormat(sample26_ATAC, from="seurat", to="anndata",assay = 'peaks',
                      outFile='../../data/scATAC-seq/sample26/Sample26_scATAC-seq_229k_processed.h5ad')
