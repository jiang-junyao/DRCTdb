library(ArchR)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
addArchRThreads(threads = 10) 
addArchRGenome("hg38")

inputFiles <- list.files('fragment',pattern = 'gz$',full.names = T)

names(inputFiles) <- map_chr(inputFiles,function(x){
  index = readLines(x,7)[7]
  if (str_starts(index,'chr')) {
    return('hg19')
  }else{
    return('hg38')
  }
})

hg38_input_file <- inputFiles[names(inputFiles) == 'hg38']


ArrowFiles <- createArrowFiles(
  inputFiles = hg38_input_file,
  sampleNames = str_extract(hg38_input_file,'(?<=t/).*(?=_f)'),
  minTSS = 4, 
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, 
  knnMethod = "UMAP",
  LSIMethod = 1
)


##Create ArchR object-----
sample27_atac <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "sample27_atac",
  copyArrows = TRUE)
sample27_atac <- filterDoublets(sample27_atac)
saveArchRProject(ArchRProj = sample27_atac, outputDirectory = "sample27_atac", load = TRUE)
#sample27_atac <- loadArchRProject('sample27_atac')
metadata <- list.files('metadata/',full.names = T) %>% map_dfr(fread,header = T)
metadata$barcode <- paste0(
  str_extract(metadata$Cell, '.*(?=#)'),
  '_atac#',
  str_extract(metadata$Cell, '(?<=#).*')
) 
cells <- intersect(metadata$barcode,getCellNames(sample27_atac)) 
metadata <- metadata %>% filter(barcode %in% cells)

sample27_atac <- sample27_atac[cells,]
sample27_atac <- addCellColData(sample27_atac,data =metadata$CellType,name = 'cell_type',cells = metadata$barcode)

sample27_atac <- addGroupCoverages(ArchRProj = sample27_atac,
                              groupBy = "cell_type",
                              useLabels = F,
                              minCells = 1000,
                              maxCells = 5000,
                              minReplicates = 3,
                              maxReplicates = 5,
                              sampleRatio = 0.9,
                              maxFragments = 50 * 10^8,
                              force = T)
sample27_atac <- addReproduciblePeakSet(
  ArchRProj = sample27_atac, 
  groupBy = "cell_type", 
  pathToMacs2 = '/home/kyh/miniconda3/envs/deeptools/bin/macs2',
  force =T,
)
sample27_atac <- addPeakMatrix(sample27_atac)
saveArchRProject(ArchRProj = sample27_atac, outputDirectory = "sample27_atac", load = TRUE)

peak_matrix <- getMatrixFromProject(
  ArchRProj = sample27_atac,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE
)
saveRDS(peak_matrix,file = 'Sample27_peak_matrix.Rds')

