library(ArchR)
setwd('LDSC_hg38/sample10/')
addArchRThreads(threads = 24) 
addArchRGenome("hg38")
inputFiles <- list.files('./fragment_file',pattern = '*.gz$',full.names = T)
ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = str_extract(inputFiles,'(?<=e/).*(?=.t)'),
    minTSS = 4, 
    minFrags = 1000, 
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
)

sample10 <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = "sample10_archR",
    copyArrows = TRUE)
sample10 <- loadArchRProject('sample10_archR/')
sample10 <- addIterativeLSI(
    ArchRProj = sample10,
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
sample10 <- addHarmony(
    ArchRProj = sample10,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force =T
)
sample10 <- addClusters(
    input = sample10,
    reducedDims = "IterativeLSI",
    method = "Seurat", 
    name = "Clusters",
    resolution = 0.8,
    force =T
)
sample10 <- addUMAP(
    ArchRProj = sample10, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force =T
)






sample10 <- addGroupCoverages(ArchRProj = sample10,
                             groupBy = "Sample",
                             minCells = 300,
                             maxCells = 30000,
                             force = TRUE)
sample10 <- addReproduciblePeakSet(
    ArchRProj = sample10, 
    groupBy = "Sample", 
    pathToMacs2 = '/home/kyh/miniconda3/bin/macs2',
    force =T,
)
sample10 <- addPeakMatrix(sample10)

saveArchRProject(ArchRProj = sample10, outputDirectory = "sample10_archR", load = TRUE)

# for i in *.gz
# do
# echo $i >> sample_info.txt
# echo zcat $i | | awk 'NR==2' >> sample_info.txt
# done

metadata <- data.table::fread('sample10_fetal_lung_metadata.txt.gz')
metadata$barcode <- str_extract(metadata$V1,'(?<=#).*')

metadata$barcode2 <- paste0(metadata$Sample,'#',metadata$barcode)
cell_barcode <- intersect(rownames(sample10@cellColData),metadata$V1) 

sample10 <- sample10[cell_barcode,]
metadata <- metadata[which(metadata$V1 %in% cell_barcode),]
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$V1
metadata <- metadata[,-1]
metadata <- metadata[rownames(sample10@cellColData),]
identical(rownames(sample10@cellColData),rownames(metadata))
sample10$cell_type <- metadata$cell_type

p1 <- plotEmbedding(ArchRProj = sample10, 
                    colorBy = "cellColData", 
                    name = "Sample", 
                    embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = sample10, 
                    colorBy = "cellColData", 
                    name = "cell_type", 
                    embedding = "UMAP")
plotPDF(p1,p2, name = "Plot-UMAP2Harmony-sample10-Clusters.pdf", ArchRProj = sample10, addDOC = FALSE, width = 5, height = 5)


saveArchRProject(ArchRProj = sample10, outputDirectory = "sample10_archR", load = TRUE)



peak_matrix <- getMatrixFromProject(
    ArchRProj = sample10,
    useMatrix = "PeakMatrix",
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE
)
saveRDS(peak_matrix,file = 'sample10_peak_matrix.Rds')
