library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
setwd('../../data/scATAC-seq/Sample8/')
addArchRThreads(threads = 1) 
addArchRGenome("hg38")
inputFiles <- list.files('./fragment',pattern = '*.gz$',full.names = T)
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = str_extract(inputFiles,'(?<=t/).*(?=_f)'),
  minTSS = 4, 
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
sample8_metadata <- data.table::fread('sample8_ATAC_metadata.txt')
sample8 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "sample8_archR",
  copyArrows = TRUE)
sample8 <- addIterativeLSI(
  ArchRProj = sample8,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( 
    resolution = c(0.2), 
    sampleCells = 30000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force =T
)
sample8 <- addHarmony(
  ArchRProj = sample8,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force =T
)
sample8 <- addClusters(
  input = sample8,
  reducedDims = "IterativeLSI",
  method = "Seurat", 
  name = "Clusters",
  resolution = 0.8,
  force =T
)
sample8 <- addUMAP(
  ArchRProj = sample8, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force =T
)
sample8 <- loadArchRProject('sample8_archR')

metadata <- read.csv('atac_barcodes.csv')

get_celltype <- function(x){
  res <- metadata[str_which(metadata$barcode,x),]$celltype
  if (length(res) > 1) {
    res <- unique(res)
    if (length(res) > 1) {
      res <- 'Unknown'
    }
  }else if(length(res) ==0 ){
    res <- 'Unknown'
  }
  return(res)
}

cell_type <- purrr::map_chr(str_extract(sample8$cellNames,'(?<=#)[A-Z]+'),get_celltype,.progress = T)
sample8 <- addCellColData(sample8,data = cell_type,name = 'cell_type',cells = sample8$cellNames)

sample8 <- sample8[sample8$cell_type != 'Unknown',]

sample8 <- addGroupCoverages(ArchRProj = sample8,
                             groupBy = "cell_type",
                             minCells = 300,
                             maxCells = 30000,
                             force = TRUE)
sample8 <- addReproduciblePeakSet(
  ArchRProj = sample8, 
  groupBy = "cell_type", 
  pathToMacs2 = '/home/kyh/miniconda3/envs/deeptools/bin/macs2',
  force =T,
)
sample8 <- addPeakMatrix(sample8)
saveArchRProject(ArchRProj = sample8, outputDirectory = "sample8_archR", load = TRUE)


peak_matrix <- getMatrixFromProject(
  ArchRProj = sample8,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE
)
saveRDS(peak_matrix,file = 'sample8_peak_matrix.Rds')