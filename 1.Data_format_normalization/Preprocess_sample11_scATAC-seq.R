library(ArchR)
setwd('.')
addArchRThreads(threads = 1) 
addArchRGenome("hg19")
inputFiles <- list.files('./fragment',pattern = '*.gz$',full.names = T)
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = paste0(rep('Control',11),1:11),
  minTSS = 4, 
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
sample11 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "sample11_archR",
  copyArrows = FALSE)
sample11 <- addIterativeLSI(
  ArchRProj = sample11,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( 
    resolution = c(0.2), 
    sampleCells = 30000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)
sample11 <- addHarmony(
  ArchRProj = sample11,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)
sample11 <- addClusters(
  input = sample11,
  reducedDims = "IterativeLSI",
  method = "Seurat", 
  name = "Clusters",
  resolution = 0.8,
  force = T
)
sample11 <- addUMAP(
  ArchRProj = sample11, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)


sample11 <- addGroupCoverages(ArchRProj = sample11,
                             groupBy = "Sample",
                             minCells = 100,
                             maxCells = 3000,
                             force = T)
sample11 <- addReproduciblePeakSet(
  ArchRProj = sample11, 
  groupBy = "Sample", 
  pathToMacs2 = '/home/kyh/miniconda3/envs/deeptools/bin/macs2',
  force =T,
)
sample11 <- addPeakMatrix(sample11)
saveArchRProject(ArchRProj = sample11, outputDirectory = "sample11_archR", load = TRUE)


peak_matrix <- getMatrixFromProject(
  ArchRProj = sample11,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE
)
saveRDS(peak_matrix,file = 'Sample11_peak_matrix.Rds')

