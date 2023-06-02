library(ArchR)
setwd('sample8/')
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

sample8 <- addGroupCoverages(ArchRProj = sample8,
                              groupBy = "Sample",
                              minCells = 300,
                              maxCells = 30000,
                              force = TRUE)
sample8 <- addReproduciblePeakSet(
  ArchRProj = sample8, 
  groupBy = "Sample", 
  pathToMacs2 = '/home/kyh/miniconda3/envs/deeptools/bin/macs2',
  force =T,
)
sample8 <- addPeakMatrix(sample8)
metadata <- sample8@cellColData |> as.data.frame()
saveArchRProject(ArchRProj = sample8, outputDirectory = "sample8_archR", load = TRUE)


peak_matrix <- getMatrixFromProject(
  ArchRProj = sample8,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE
)
saveRDS(peak_matrix,file = 'sample8_peak_matrix.Rds')