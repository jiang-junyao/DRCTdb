library(ArchR)

setwd('LDSC_hg38/sample3/')
addArchRThreads(threads = 40) 
addArchRGenome("hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
inputFiles <- list.files('./fragment_file',pattern = '*.gz$',full.names = T)
ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = str_extract(inputFiles,'(?<=e/).*(?=_)'),
    minTSS = 4, 
    minFrags = 1000, 
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
)
sample3 <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = "sample3_archR",
    copyArrows = TRUE)
sample3 <- loadArchRProject(path = "sample3_archR")
sample3 <- addIterativeLSI(
    ArchRProj = sample3,
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
sample3 <- addHarmony(
    ArchRProj = sample3,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force =T
)

sample3 <- addClusters(
    input = sample3,
    reducedDims = "IterativeLSI",
    method = "Seurat", 
    name = "Clusters",
    resolution = 0.8,
    force =T
)
sample3 <- addUMAP(
    ArchRProj = sample3, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force =T
)
sample3 <- addGroupCoverages(ArchRProj = sample3,
                              groupBy = "Sample",
                              minCells = 300,
                              maxCells = 30000,
                              force = TRUE)
sample3 <- addReproduciblePeakSet(
    ArchRProj = sample3, 
    groupBy = "Sample", 
    pathToMacs2 = '/home/kyh/miniconda3/bin/macs2',
    force =T,
)
sample3 <- addPeakMatrix(sample3)
saveArchRProject(ArchRProj = sample3, outputDirectory = "sample3_archR", load = TRUE)
peak_matrix <- getMatrixFromProject(
    ArchRProj = sample3,
    useMatrix = "PeakMatrix",
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE
)
saveRDS(peak_matrix,file = 'sample3_peak_matrix.Rds')
