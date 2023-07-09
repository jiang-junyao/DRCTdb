library(ArchR)
library(purrr)
options(timeout = max(1000, getOption("timeout"))) 
addArchRGenome("hg38")
addArchRThreads(threads = 1) 
saveplot <- function(object,filenames = filename,width = NULL,height = NULL,dpi = 600){
  if (is.null(width)| is.null(height)) {
    stop('Must input figure width and height')
  }
  
  if(!dir.exists("Fig")){
    dir.create("Fig")
  }
  map(c("TIFF","PNG","PDF"),function(x){
    if(!dir.exists(paste0('Fig/',x))){
      dir.create(paste0('Fig/',x))
    }
  })
  
  pdf_filenames = paste0('Fig/PDF/',filenames,'.pdf')
  pdf_width = width*2
  pdf_height = height*2
  pdf_res = dpi
  
  tiff_filenames = paste0('Fig/TIFF/',filenames,'.tiff')
  tiff_width = width*1000
  tiff_height = height*1000
  tiff_res = dpi
  
  png_filenames = paste0('Fig/PNG/',filenames,'.png')
  png_width = width*500
  png_height = height*500
  if (dpi > 600) {
    png_res = 300
  }else{
    png_res = dpi/2
  }
  
  
  cairo_pdf(file = pdf_filenames,width = pdf_width, height = pdf_height,fallback_resolution = pdf_res)
  plot(object)
  dev.off()
  
  tiff(filename = tiff_filenames,width = tiff_width, height = tiff_height, units = "px", res = tiff_res, compression = "lzw")
  plot(object)
  dev.off()
  
  png(filename = png_filenames,width = png_width, height = png_height, res = png_res,units = "px")
  plot(object)
  dev.off()
  
}



##-----
ArrowFiles <- createArrowFiles(
  inputFiles = 'GSE214979_atac_fragments.tsv.gz',
  sampleNames = 'Sample24',
  minTSS = 3, 
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
sample24 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "sample24_archR",
  copyArrows = TRUE)

metadata <- data.table::fread('GSE214979_cell_metadata.csv.gz')
metadata$barcode <- paste0('Sample24#',metadata$V1)
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$barcode
cell_names <- intersect(metadata$barcode,getCellNames(sample24))
metadata <- metadata[cell_names,]

sample24 <- sample24[cell_names,]

identical(metadata$barcode,getCellNames(sample24))

for (i in colnames(metadata)[-1]) {
  sample24 <- addCellColData(sample24,data = metadata[[i]],name = i,cells = rownames(metadata))
}

sample24$cell_type <- sample24$subs


sample24_rna_h5 <- import10xFeatureMatrix(
  input = 'GSE214979_filtered_feature_bc_matrix.h5',
  names = 'Sample24'
)
cell_names <- intersect(colnames(sample24_rna_h5),getCellNames(sample24)) 
sample24_rna_h5 <- sample24_rna_h5[,cell_names]
sample24 <- sample24[cell_names,]
identical(colnames(sample24_rna_h5),getCellNames(sample24))

sample24 <- addGeneExpressionMatrix(input = sample24,
                                    seRNA = sample24_rna_h5,
                                    force = TRUE)
saveRDS(sample24_rna_h5,file = 'Sample24_exp_matrix.Rds')

sample24 <- addIterativeLSI(
  ArchRProj = sample24, 
  clusterParams = list(
    resolution = 0.3, 
    sampleCells = 30000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  sampleCellsPre = 60000,
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA"
)
sample24 <- addHarmony(
  ArchRProj = sample24,
  reducedDims = "LSI_RNA",
  name = "Harmony",
  groupBy = "id",
  corCutOff = 0.8,
  force = TRUE
)
sample24 <- addUMAP(
  ArchRProj = sample24, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony_RNA", 
  nNeighbors = 30, 
  minDist = 0.3, 
  metric = "cosine",
  force = TRUE
)
p1 <- plotEmbedding(ArchRProj = sample24, colorBy = "cellColData", name = "cell_type", embedding = "UMAPHarmony_RNA")
saveplot(p1,filenames = 'Plot-UMAP2Harmony-cell_type_RNA',width = 5,height = 5)
plotPDF(p1, name = "Plot-UMAP2Harmony-cell_type_RNA.pdf", ArchRProj = sample24, addDOC = FALSE, width = 5, height = 5)

sample24 <- addGroupCoverages(ArchRProj = sample24,
                              groupBy = "cell_type",
                              minCells = 20000,
                              maxCells = 50000,
                              maxReplicates = 5,
                              force = T)
sample24 <- addReproduciblePeakSet(
  ArchRProj = sample24, 
  groupBy = "cell_type", 
  pathToMacs2 = '/home/kyh/miniconda3/envs/RNA-seq/envs/deeptools/bin/macs2',
  force =T
)
sample24 <- addPeakMatrix(sample24)
saveArchRProject(ArchRProj = sample24, outputDirectory = "sample24_archR", load = TRUE)

sample24_peak_matrix <- getMatrixFromProject(
  ArchRProj = sample24,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE
)
saveRDS(sample24_peak_matrix,file = 'Sample24_peak_matrix.Rds')
##
library(Signac)
library(Seurat)
library(tidyverse)
source('../2.Data preprocessing/preprocess_functions.R')
sample24_ATAC <- readRDS('../../data/scATAC-seq/Sample24/Sample24_peak_matrix.Rds')


coldata <- as.data.frame(sample24_ATAC@colData)
coldata$raw_barcode <- rownames(coldata)
sparse_mtx <- sample24_ATAC@assays@data$PeakMatrix
rownames(sparse_mtx) <- paste(as.data.frame(sample24_ATAC@rowRanges)[[1]],as.data.frame(sample24_ATAC@rowRanges)[[2]],as.data.frame(sample24_ATAC@rowRanges)[[3]],sep = '-')
sparse_mtx <- sparse_mtx[which(map_vec(rownames(sparse_mtx),subset_peaks)),]

chrom_assay <- CreateChromatinAssay(
  counts = sparse_mtx,
  sep = c("-", "-"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)
sample24 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

all(colnames(sample24) %in% colnames(sample24_ATAC))
identical(colnames(sample24),colnames(sample24_ATAC))

metadata <- sample24_ATAC@colData |> as.data.frame()
identical(colnames(sample24),rownames(metadata))

sample24 <- AddMetaData(sample24,metadata = metadata)
saveRDS(sample24,file = '../../data/scATAC-seq/Sample24/sample24_scATAC-seq_100k_all_processed.Rds')

sample24_healthy <- subset(x = sample24,Diagnosis=='Unaffected')

saveRDS(sample24_healthy,file = '../../data/scATAC-seq/Sample24/sample24_scATAC-seq_100k_healthy_processed.Rds')



