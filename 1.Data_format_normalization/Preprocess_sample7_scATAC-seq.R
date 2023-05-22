library(ArchR)
setwd('sample7/')
addArchRThreads(threads = 1) 
addArchRGenome("hg38")

inputFiles <- list.files('./fragment_file',pattern = '*.gz$',full.names = T)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = str_extract(inputFiles,'(?<=73_).*(?=_)'),
  minTSS = 4, 
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

sample7 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "sample7_archR",
  copyArrows = TRUE)
#sample7 <- loadArchRProject('sample7_archR/',showLogo = F)
sample7 <- addIterativeLSI(
  ArchRProj = sample7,
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
sample7 <- addHarmony(
  ArchRProj = sample7,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)
sample7 <- addClusters(
  input = sample7,
  reducedDims = "IterativeLSI",
  method = "Seurat", 
  name = "Clusters",
  resolution = 0.8
)
sample7 <- addUMAP(
  ArchRProj = sample7, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = sample7, 
                    colorBy = "cellColData", 
                    name = "Sample", 
                    embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = sample7, 
                    colorBy = "cellColData", 
                    name = "Clusters", 
                    embedding = "UMAP")
plotPDF(p1,p2, name = "Plot-UMAP2Harmony-Sample7-Clusters.pdf", ArchRProj = sample7, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = sample7, outputDirectory = "sample7_archR", load = FALSE)

sample7 <- addGroupCoverages(ArchRProj = sample7,
                             groupBy = "Sample",
                             minCells = 300,
                             maxCells = 10000,
                             force = T)
sample7 <- addReproduciblePeakSet(
  ArchRProj = sample7, 
  groupBy = "Sample", 
  pathToMacs2 = '/home/kyh/miniconda3/envs/deeptools/bin/macs2',
  force =T,
)
sample7 <- addPeakMatrix(sample7)
saveArchRProject(ArchRProj = sample7, outputDirectory = "sample7_archR", load = TRUE)


peak_matrix <- getMatrixFromProject(
  ArchRProj = sample7,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE
)
saveRDS(peak_matrix,file = 'Sample7_peak_matrix.Rds')
