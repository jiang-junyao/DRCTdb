library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg19)
setwd('~/sample26/')
addArchRThreads(threads = 1) 
addArchRGenome("hg19")
# df <- data.table::fread('fragment_file/GSM3722043_BY_10k_fragments.tsv.gz')
# df <- df %>% dplyr::filter(str_detect(V1,'hg19'))
# df$V1 <- str_extract(df$V1,'chr.*')
# data.table::fwrite(df,'fragment_file/GSM3722043_BY_10k_fragments.tsv',sep = '\t',col.names = F)
# df <- data.table::fread('fragment_file/GSM3722046_BY_5k_fragments.tsv.gz')
# df <- df %>% dplyr::filter(str_detect(V1,'hg19'))
# df$V1 <- str_extract(df$V1,'chr.*')
# data.table::fwrite(df,'fragment_file/GSM3722046_BY_5k_fragments.tsv',sep = '\t',col.names = F)
# df <- data.table::fread('fragment_file/GSM3722044_BY_1k_fragments.tsv.gz')
# df <- df %>% dplyr::filter(str_detect(V1,'hg19'))
# df$V1 <- str_extract(df$V1,'chr.*')
# data.table::fwrite(df,'fragment_file/GSM3722044_BY_1k_fragments.tsv',sep = '\t',col.names = F)
# df <- data.table::fread('fragment_file/GSM3722045_BY_500_fragments.tsv.gz')
# df <- df %>% dplyr::filter(str_detect(V1,'hg19'))
# df$V1 <- str_extract(df$V1,'chr.*')
# data.table::fwrite(df,'fragment_file/GSM3722045_BY_500_fragments.tsv',sep = '\t',col.names = F)

Optimization <- function(iterations,sampleCells,varFeatures,archr_obj = NULL,groub_by = 'Sample'){
  archr_obj <- addIterativeLSI(
    ArchRProj = archr_obj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = iterations, 
    clusterParams = list( 
      resolution = c(0.1, 0.8,1.5, 2)[1:(iterations-1)], 
      sampleCells = sampleCells, 
      n.start = 10
    ), 
    sampleCellsPre = sampleCells,
    varFeatures = varFeatures, 
    saveIterations = F,
    dimsToUse = 1:50,
    force = TRUE 
  )
  
  archr_obj <- addHarmony(
    ArchRProj = archr_obj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    corCutOff = 0.8,
    force = TRUE
  )
  
  archr_obj <- addUMAP(
    ArchRProj = archr_obj, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.3, 
    metric = "cosine",
    force = TRUE
  )
  return(archr_obj)
}
inputFiles <- list.files('./fragment_file',pattern = '*.gz$',full.names = T)
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = str_extract(inputFiles,'(?<=e/).*(?=.t)'),
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
ArrowFiles <- list.files('.',pattern = '*.arrow',full.names = T)
sample26 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Sample26_archR",
    copyArrows = TRUE)

sample26 <- sample26 <- filterDoublets(sample26)
saveArchRProject(sample26,outputDirectory = 'Sample26_archR/')
sample26 <- Optimization(iterations = 5,sampleCells = 10000,varFeatures = 50000,archr_obj = sample26,groub_by = "Sample")
saveArchRProject(ArchRProj = sample26, outputDirectory = "Sample26_archR", load = TRUE)
p1 <- plotEmbedding(ArchRProj = sample26, colorBy = "cellColData", name = 'Sample', embedding = "UMAPHarmony")
sample26 <- loadArchRProject('Sample26_archR')
#load metadata
meta_files <- list.files('./metadata/',full.names = T)
metadatas  <- map_dfr(
  meta_files,function(x){
    data.table::fread(x,header = T)
  }
)
metadatas <- metadatas %>% group_by(Group_Barcode) %>% summarise(across(everything(),~ head(.x, n =1)))


coldata <- getCellColData(sample26) |> as.data.frame()
coldata$cell_name <- rownames(coldata) 
coldata$barcode <- str_extract(rownames(coldata),'(?<=_).*') |> str_replace('_fragments','')
coldata <- coldata %>% left_join(metadatas,by = c('barcode' = 'Group_Barcode'))
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$cell_name 
coldata <- coldata[getCellNames(sample26),]
identical(rownames(coldata),getCellNames(sample26))

sample26 <- addCellColData(sample26,data = coldata$Clusters,name = 'annotated_cluster',cells = rownames(coldata))
sample26 <- addCellColData(sample26,data = coldata$Group,name = 'Group',cells = rownames(coldata))
sample26 <- addCellColData(sample26,data = coldata$Internal_Name,name = 'Internal_Name',cells = rownames(coldata))
saveArchRProject(ArchRProj = sample26, outputDirectory = "Sample26_archR", load = TRUE)




library(AtacAnnoR)
sample26 <- loadArchRProject('Sample26_archR')
metadata <- as.data.frame(sample26@cellColData)

metadata_1 <- metadata[!is.na(metadata$Group),]
metadata_2 <- metadata[is.na(metadata$Group),]
metadata_2$Group <- str_extract(metadata_2$Sample,'(?<=\\d_).*(?=_fragment)')
metadata <- rbind(metadata_1,metadata_2)

sample26 <- addCellColData(sample26,cells = rownames(metadata),data = metadata$Group,name = 'Group',force = T)

#annotate pbmc cells
pbmc_data <- sample26[which(str_detect(sample26$Group,'PBMC|BY|pbmc')),]
pbmc_ref <- readRDS('pbmc_multimodal_2023.rds')
pbmc_ref <- subset(pbmc_ref,cells = names(which(pbmc_ref$celltype.l2 != 'Doublet')))
pbmc_data <- RunAtacAnnoR_ArchR(
  query_ArchRproj = pbmc_data,
  ref_mtx = pbmc_ref[['SCT']]@counts,
  ref_celltype = pbmc_ref$celltype.l2
)
  
pbmc_anno <- pbmc_data$final_pred
names(pbmc_anno) <- pbmc_data$cellNames
#annotate bone marrow cells

bone_data <- sample26[which(str_detect(sample26$Group,'arrow')),]
bone_ref <- readRDS('sample16_Bone_marrow_RNA_Healthy_35k_processed.Rds')
bone_data <- RunAtacAnnoR_ArchR(
    query_ArchRproj = bone_data,
    ref_mtx = bone_ref[['RNA']]@counts,
    ref_celltype = bone_ref$cell_type
)

bone_anno <- bone_data$final_pred
names(bone_anno) <- bone_data$cellNames


temp_cells <- sample26$Group[which(!(sample26$cellNames %in% c(names(pbmc_anno),names(bone_data))))]
names(temp_cells) <- sample26$cellNames[which(!(sample26$cellNames %in% c(names(pbmc_anno),names(bone_data))))]
get_celltype <- function(x){
    if (x == '0p1_99p9_CD4Mem_CD8Naive') {
        return('CD4Mem CD8Naive')
    }else if(x == '0p1_99p9_Mono_T'){
        return('Mono T')
    }else if(x == '0p5_99p5_CD4Mem_CD8Naive'){
        return('CD4Mem CD8Naive')
    }else if(x == '0p5_99p5_Mono_T'){
        return('Mono T')
    }else if(x == '1_99_CD4Mem_CD8Naive'){
        return('CD4Mem CD8Naive')
    }else if(x == '1_99_Mono_T'){
        return('Mono T')
    }else if(x == '50_50_CD4Mem_CD8Naive'){
        return('CD4Mem CD8Naive')
    }else if(x == '50_50_Mono_T'){
        return('Mono T')
    }else if(x == '99_1_CD4Mem_CD8Naive'){
        return('CD4Mem CD8Naive')
    }else if(x == '99_1_MonoT'){
        return('Mono T')
    }else if(x == '99p9_0p1_Mono_T'){
        return('Mono T')
    }else if(x == '99p5_0p5_MonoT'){
        return('Mono T')
    }else if(x == '99p5_0p5_CD4Mem_CD8Naive'){
        return('CD4Mem CD8Naive')
    }else if(x == '99p9_0p1_CD4Mem_CD8Naive'){
        return('Mono T')
    }else if(x == '99p5_0p5_CD4Mem_CD8Naive'){
        return('Mono T')
    }else{
        return(x)
    }
}
new_celltemp <- purrr::map_chr(temp_cells,get_celltype)

cell_type <- c(new_celltemp,pbmc_anno,bone_anno)
sample26 <- addCellColData(sample26,cells = names(cell_type),data = cell_type,name = 'cell_type')
saveArchRProject(ArchRProj = sample26, outputDirectory = "Sample26_archR", load = TRUE)
sample26 <- loadArchRProject('Sample26_archR')
sample26 <- addGroupCoverages(
    ArchRProj = sample26,
    groupBy = "cell_type",
    minCells = 1000,
    maxCells = 5000,
    minReplicates = 3,
    maxReplicates = 5,
    force = T
)
sample26 <- addReproduciblePeakSet(
    ArchRProj = sample26, 
    groupBy = "cell_type", 
    peaksPerCell = 2000,
    sampleLabels = NULL,
    cutOff = 0.05,
    method = "q",
    maxPeaks = 400000,
    pathToMacs2 = '/home/kyh/miniconda3/envs/deeptools/bin/macs2',
    force =T
)
sample26 <- addPeakMatrix(sample26)
saveArchRProject(ArchRProj = sample26, outputDirectory = "Sample26_archR", load = TRUE)

peak_matrix <- getMatrixFromProject(
  ArchRProj = sample26,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE
)
saveRDS(peak_matrix,file = 'Sample26_peak_matrix.Rds')
