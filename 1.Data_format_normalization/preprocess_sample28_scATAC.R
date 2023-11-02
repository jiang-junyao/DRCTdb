library(ArchR)
library(dplyr)
library(tidyr)
library(mclust)
library(ggrastr)
##function
identifyCells <- function(df, TSS_cutoff=6, nFrags_cutoff=2000, minTSS=5, minFrags=1000, maxG=4){
  # Identify likely cells based on gaussian mixture modelling.
  # Assumes that cells, chromatin debris, and other contaminants are derived from
  # distinct gaussians in the TSS x log10 nFrags space. Fit a mixture model to each sample
  # and retain only cells that are derived from a population with mean TSS and nFrags passing
  # cutoffs
  ####################################################################
  # df = data.frame of a single sample with columns of log10nFrags and TSSEnrichment
  # TSS_cutoff = the TSS cutoff that the mean of a generating gaussian must exceed
  # nFrags_cutoff = the log10nFrags cutoff that the mean of a generating gaussian must exceed
  # minTSS = a hard cutoff of minimum TSS for keeping cells, regardless of their generating gaussian
  # maxG = maximum number of generating gaussians allowed
  
  cellLabel <- "cell"
  notCellLabel <- "not_cell"
  
  if(nFrags_cutoff > 100){
    nFrags_cutoff <- log10(nFrags_cutoff)
    minFrags <- log10(minFrags)
  } 
  
  # Fit model
  set.seed(1)
  mod <- Mclust(df, G=2:maxG, modelNames="VVV")
  
  # Identify classifications that are likely cells
  means <- mod$parameters$mean
  
  # Identify the gaussian with the maximum TSS cutoff
  idents <- rep(notCellLabel, ncol(means))
  idents[which.max(means["TSSEnrichment",])] <- cellLabel
  
  names(idents) <- 1:ncol(means)
  
  # Now return classifications and uncertainties
  df$classification <- idents[mod$classification]
  df$classification[df$TSSEnrichment < minTSS] <- notCellLabel
  df$classification[df$nFrags < minFrags] <- notCellLabel
  df$cell_uncertainty <- NA
  df$cell_uncertainty[df$classification == cellLabel] <- mod$uncertainty[df$classification == cellLabel]
  return(list(results=df, model=mod))
}
# Helper functions for ArchR


getMatrixValuesFromProj <- function(proj, matrixName="GeneScoreMatrix", names=NULL, imputeMatrix=FALSE){
  # Return the imputed matrix from an ArchR project
  # Must have already added imputedWeights, etc.
  # Names is a vector of feature names to return. If not provided, will use all available features
  # Warning though: imputing the matrix for all features may take a very long time
  
  # Returns a summarized experiment:
  se <- getMatrixFromProject(proj, useMatrix = matrixName, binarize = FALSE)
  
  # Get mat with cell names and row names
  mat <- assays(se)[[matrixName]]
  colnames(mat) <- rownames(colData(se))
  rownames(mat) <- rowData(se)$name # All matrix rowData has a name field
  
  # Subset by provided names
  if(!is.null(names)){
    # Check if any names are invalid
    validNames <- names[names %in% rownames(mat)]
    if(any(!names %in% validNames)){
      invalidNames <- names[!names %in% validNames]
      message(sprintf("Warning! name(s) %s are not present in matrix!", paste(invalidNames, collapse=',')))
    }
    mat <- mat[validNames,]
  }
  
  # Impute matrix values
  if(imputeMatrix){
    message("Imputing matrix...")
    imputeWeights <- getImputeWeights(proj)
    mat <- ArchR::imputeMatrix(mat = as.matrix(mat), imputeWeights = imputeWeights)
  }
  mat
}


getClusterPeaks <- function(proj, clusterNames, peakGR=NULL, replicateScoreQuantileCutoff=0, originalScore=FALSE){
  # This function will return the subset of peaks from the full ArchR project that were 
  # initially called on the clusters provided in clusterNames.
  ######################################################################################
  # proj = ArchR project
  # clusterNames = name or names of clusters to pull peaks from. These cluster names must
  #   match the cluster names originally used to call peaks
  # peakGR = the ArchR peak genomic range obtained using 'getPeakSet'. Will return a subset of
  #   these peaks that overlap the peaks originally called using the provided clusters
  # replicateScoreQuantileCutoff = A numeric quantile cutoff for selecting peaks. 
  if(is.null(peakGR)){
    peakGR <- getPeakSet(proj)
  }
  peakDir <- paste0(proj@projectMetadata$outputDirectory, "/PeakCalls")
  calledPeaks <- lapply(clusterNames, function(x){
    readRDS(paste0(peakDir, sprintf("/%s-reproduciblePeaks.gr.rds", x)))
  }) %>% as(., "GRangesList") %>% unlist()
  calledPeaks <- calledPeaks[calledPeaks$replicateScoreQuantile >= replicateScoreQuantileCutoff]
  peakGR <- peakGR[overlapsAny(peakGR, calledPeaks)]
  if(originalScore){
    message(sprintf("Getting original scores from clusters..."))
    # Replace 'score' column with the score of the original peak call in this cluster
    # (if multiple clusters, replaces with the maximum score)
    ol <- findOverlaps(peakGR, calledPeaks, type="any", maxgap=0, ignore.strand=TRUE)
    odf <- as.data.frame(ol)
    odf$og_score <- calledPeaks$score[odf$subjectHits]
    score_df <- odf %>% group_by(queryHits) %>% summarize(max_score=max(og_score)) %>% as.data.frame()
    peakGR$score[score_df$queryHits] <- score_df$max_score
  }
  peakGR
}


buildUMAPdfFromArchR <- function(proj, cellColData=NULL, embeddingName="UMAP", 
                                 useCells=NULL, dataMat=NULL, featureName=NULL, shuffle=TRUE, 
                                 lowerPctLim=NULL, upperPctLim=NULL){
  # Return a three column UMAP df from an ArchR project
  # If cellColData is not null, return the indicated column
  # dataMat is a pre-populated cell x feature matrix of values to plot. 
  # The featureName indicates which one
  
  # Get UMAP coordinates first:
  df <- proj@embeddings[[embeddingName]]$df
  if(is.null(useCells)){
    useCells <- rownames(df)
  }
  colnames(df) <- c("UMAP1", "UMAP2")
  df <- df[useCells,] %>% as.data.frame()
  if(!is.null(cellColData)){
    df[,3] <- proj@cellColData[useCells,cellColData] %>% as.vector()
    colnames(df) <- c("UMAP1", "UMAP2", cellColData)
  }
  if(!is.null(dataMat) & !is.null(featureName)){
    df <- merge(df, dataMat, by=0, all=TRUE)
    df <- df[,c("UMAP1", "UMAP2", featureName)] 
  }
  if(shuffle){
    df <- df[sample(nrow(df), replace=FALSE),]
  }
  # Force limits if indicated
  if(!is.null(lowerPctLim)){
    lowerLim <- quantile(df[,3], probs=c(lowerPctLim))
    df[,3][df[,3] <= lowerLim] <- lowerLim
  }
  if(!is.null(upperPctLim)){
    upperLim <- quantile(df[,3], probs=c(upperPctLim))
    df[,3][df[,3] >= upperLim] <- upperLim
  }
  df
}


scoreGeneSet <- function(expr, geneSet){
  # Generate scores for each cell in expr matrix (log2TP10K, genes x cells)
  # See: Smillie et al. Cell 2019
  
  # Subset expr matrix by genes in geneSet:
  validGenes <- geneSet[geneSet %in% rownames(expr)]
  subExpr <- expr[validGenes,]
  
  # Remove any genes that have no expression in any cells
  subExpr <- subExpr[rowSums(subExpr) > 0,]
  
  # Prevent highly expressed genes from dominating gene score signature by
  # scaling each gene by its root mean squared expression
  scaledSubExpr <- subExpr %>% t() %>% scale(., center=FALSE) %>% t()
  
  # Signature score is the mean scaled expression across all genes in signature
  scores <- colMeans(scaledSubExpr)
  return(scores)
}


# Functions for creating 'low-overlapping aggregates' of cells

computeKNN <- function(data=NULL, query=NULL, k=50, includeSelf=FALSE, ...){
  # Compute KNN for query points (usually a reduced dims matrix)
  # This returns a matrix of indices mapping query to neighbors in data
  # If query has n cells (rows) and k = 50, will be a n x 50 matrix
  if(is.null(query)){
    query <- data
    searchSelf <- TRUE
  }else{
    searchSelf <- FALSE
  }
  if(searchSelf & !includeSelf){
    knnIdx <- nabor::knn(data = data, query = query, k = k + 1, ...)$nn.idx
    knnIdx <- knnIdx[,-1,drop=FALSE]
  }else{
    knnIdx <- nabor::knn(data = data, query = query, k = k, ...)$nn.idx
  }
  knnIdx
}


getLowOverlapAggregates <- function(proj, target.agg=500, k=100, overlapCutoff=0.8, dimReduc="IterativeLSI", seed=1){
  # Generate low-overlapping aggregates of cells
  ##############################################
  # proj = ArchR project
  # target.agg = number of target aggregates (before filtering)
  # k = number of cells per aggreagate
  # overlapCutoff = Maximum allowable overlap between aggregates
  set.seed(seed)
  
  # Get reduced dims:
  rD <- getReducedDims(proj, reducedDims=dimReduc)
  
  # Subsample
  idx <- sample(seq_len(nrow(rD)), target.agg, replace = !nrow(rD) >= target.agg)
  
  # Get KNN Matrix:
  knnObj <- computeKNN(data=rD, query=rD[idx,], k=k)
  
  # Check whether aggregates pass overlap cutoff
  keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * k))
  
  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]
  
  # Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList
  
  # Name aggregates and return as a df of cell ids x aggs
  names(knnObj) <- paste0("agg", seq_len(length(knnObj)))
  knnDF <- data.frame(knnObj)[,c(3,2)]
  colnames(knnDF) <- c("cell_name", "group")
  knnDF$cell_name <- as.character(knnDF$cell_name)
  knnDF
}


# Cluster visualization helpers

relabelClusters <- function(proj, clusterName="Clusters"){
  # Relabel clusters to be ordered by cluster size
  
  ogClusters <- getCellColData(proj)[[clusterName]]
  tabDF <- base::table(ogClusters) %>% as.data.frame
  colnames(tabDF) <- c("Clusters", "count")
  tabDF["NewClusters"] <- rank(-tabDF$count)
  swapVec <- paste0("C", tabDF$NewClusters)
  names(swapVec) <- tabDF$Clusters
  
  # Now replace cluster names
  newClust <- sapply(ogClusters, function(x) swapVec[x]) %>% unname()
  proj <- addCellColData(proj, data=newClust, name=clusterName, cells=getCellNames(proj), force=TRUE)
  return(proj)
}


visualizeClustering <- function(proj, pointSize=0.75, prefix="", clusterName="Clusters", sampleName="Sample2", embedding="UMAP", 
                                sampleCmap=NULL, diseaseCmap=NULL, barwidth=0.9){
  # Plot various clustering results
  # Set colormap
  qualcmap <- cmaps_BOR$stallion
  quantcmap <- cmaps_BOR$solarExtra
  namedSampCmap <- TRUE
  namedDiseaseCmap <- TRUE
  
  if(is.null(sampleCmap)){
    sampleCmap <- qualcmap
    namedSampCmap <- FALSE
  }
  
  if(is.null(diseaseCmap)){
    diseaseCmap <- qualcmap
    namedDiseaseCmap <- FALSE
  }
  
  # Plot the UMAPs by Sample and Cluster:
  p1 <- plotEmbedding(proj, colorBy="cellColData", name=sampleName, embedding=embedding, plotAs="points", size=pointSize, pal=sampleCmap, labelMeans=FALSE)
  p2 <- plotEmbedding(proj, colorBy="cellColData", name=clusterName, embedding= embedding, plotAs="points", size=pointSize, labelMeans=FALSE)
  p3 <- plotEmbedding(proj, colorBy="cellColData", name="diseaseStatus", embedding=embedding, plotAs="points", size=pointSize, pal=diseaseCmap, labelMeans=FALSE)
  proj@cellColData$log10nFrags <- log10(proj@cellColData$nFrags)
  p4 <- plotEmbedding(proj, colorBy="cellColData", name="log10nFrags", embedding=embedding, plotAs="points", size=pointSize, labelMeans=FALSE)
  p5 <- plotEmbedding(proj, colorBy="cellColData", name="TSSEnrichment", embedding=embedding, plotAs="points", size=pointSize, labelMeans=FALSE)
  p6 <- plotEmbedding(proj, colorBy="cellColData", name="DoubletScore", 
                      embedding = embedding, plotAs="points", size=pointSize, labelMeans=FALSE, imputeWeights=getImputeWeights(proj))
  p7 <- plotEmbedding(proj, colorBy = "cellColData", name="cellCallUncertainty", 
                      embedding = embedding, plotAs="points", size=pointSize, labelMeans=FALSE, imputeWeights=getImputeWeights(proj))
  ggAlignPlots(p1,p2,p3,p4,p5,p6,p7, type="h")
  plotPDF(p1,p2,p3,p4,p5,p6,p7, name = paste0(prefix,"Plot-UMAP-Sample-Clusters.pdf"), ArchRProj=proj, addDOC=FALSE, width=5, height=5)
  
  # Non-ArchR plots:
  plotDir <- paste0(proj@projectMetadata$outputDirectory, "/Plots")
  
  # Bar plot cluster counts
  clustVec <- getCellColData(proj)[[clusterName]] %>% gsub("[^[:digit:].]", "", .) %>% as.numeric()
  tabDF <- base::table(clustVec) %>% as.data.frame
  colnames(tabDF) <- c("Clusters", "count")
  
  pdf(paste0(plotDir,sprintf("/%sclusterBarPlot.pdf", prefix)))
  print(qcBarPlot(tabDF, cmap=qualcmap, barwidth=barwidth))
  dev.off()
  
  # Stacked bar plot fraction samples in clusters
  clustBySamp <- fractionXbyY(clustVec, proj$Sample2, add_total=TRUE, xname="Cluster", yname="Sample")
  
  pdf(paste0(plotDir, sprintf("/%sclustBySampleBarPlot.pdf", prefix)))
  print(stackedBarPlot(clustBySamp, cmap=sampleCmap, namedColors=namedSampCmap, barwidth=barwidth))
  dev.off()
  
  # Stacked bar plot fraction disease in clusters
  diseaseBySamp <- fractionXbyY(clustVec, proj$diseaseStatus, add_total=TRUE, xname="Cluster", yname="diseaseStatus")
  
  pdf(paste0(plotDir, sprintf("/%sclustByDiseaseBarPlot.pdf", prefix)))
  print(stackedBarPlot(diseaseBySamp, cmap=diseaseCmap, namedColors=namedDiseaseCmap, barwidth=barwidth))
  dev.off()
  
  return(proj)
}


# Functions for working with peak to gene linkages

getP2G_GR <- function(proj, corrCutoff=NULL, varCutoffATAC=0.25, varCutoffRNA=0.25, filtNA=TRUE){
  # Function to get peaks and genes involved in peak to gene links
  # (See: https://github.com/GreenleafLab/ArchR/issues/368)
  ############################################################
  # proj: ArchR project that alreayd has Peak2GeneLinks
  # corrCutoff: minimum numeric peak-to-gene correlation to return
  # varCutoffATAC: minimum variance quantile of the ATAC peak accessibility when selecting links
  # varCutoffRNA: minimum variance quantile of the RNA gene expression when selecting links
  p2gDF <- metadata(proj@peakSet)$Peak2GeneLinks
  p2gDF$symbol <- mcols(metadata(p2gDF)$geneSet)$name[p2gDF$idxRNA] %>% as.character()
  p2gDF$peakName <- (metadata(p2gDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2gDF$idxATAC]
  # Remove peaks with 'NA' correlation values
  if(filtNA){
    p2gDF <- p2gDF[!is.na(p2gDF$Correlation),]
  }
  if(!is.null(corrCutoff)){
    p2gDF <- p2gDF[(p2gDF$Correlation > corrCutoff),]
  }
  # Filter by variance quantile
  p2gDF <- p2gDF[which(p2gDF$VarQATAC > varCutoffATAC & p2gDF$VarQRNA > varCutoffRNA),]
  # The genomic range contains just the peak ranges:
  p2gGR <- metadata(p2gDF)$peakSet[p2gDF$idxATAC]
  mcols(p2gGR) <- p2gDF
  p2gGR
}


grLims <- function(gr){
  # Get the minimum and maximum range from a GR
  if(length(gr) == 0){
    return(NA)
  }
  starts <- start(gr)
  ends <- end(gr)
  c(min(starts, ends), max(starts, ends))
}


getP2Gregions <- function(proj, genes, p2gGR=NULL, corrCutoff=0.4, buffer_space=0.05, min_width=25000, ...) {
  # Function to get regions containing entire peak to gene region,
  # i.e. a GR that contains all peak to gene links
  ###############################################################
  # p2gGR: genomic range containing all peak to gene links
  # genes: vector of genes to look up
  # buffer_space: fraction of total length to expand on each side of region
  
  # Get gene GR from ArchR project
  geneGR <- promoters(getGenes(proj)) # Promoters gets 2kb upstream and 200bp downstream
  geneGR <- geneGR[!is.na(geneGR$symbol)]
  
  # if p2gGR not provided, pull it from ArchR project
  if(is.null(p2gGR)){
    p2gGR <- getP2G_GR(proj, corrCutoff=corrCutoff, ...)
  }
  
  # Now for each gene, construct GR of all loops and gene TSS
  resultGR <- geneGR[match(genes, geneGR$symbol)]
  start(resultGR) <- sapply(resultGR$symbol, function(x){
    min(grLims(resultGR[resultGR$symbol == x]), grLims(p2gGR[p2gGR$symbol == x]), na.rm=TRUE)
  })
  end(resultGR) <- sapply(resultGR$symbol, function(x){
    max(grLims(resultGR[resultGR$symbol == x]), grLims(p2gGR[p2gGR$symbol == x]), na.rm=TRUE)
  })
  
  # Finally, resize by buffer space
  resultGR <- resize(resultGR, width=width(resultGR) + buffer_space*width(resultGR), fix="center")
  resultGR <- resize(resultGR, width=ifelse(width(resultGR) > min_width, width(resultGR), min_width), fix="center")
  resultGR
}



##------
addArchRThreads(threads = 6)
addArchRGenome('hg38')
inputFiles <- getInputFiles("fragment/")
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 0, 
  minFrags = 1000, 
  addTileMat = FALSE, 
  addGeneScoreMat = FALSE
)
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "unfiltered_output"
)


minTSS <- 5
samples <- unique(proj$Sample)
cellData <- getCellColData(proj)
cellResults <- lapply(samples, function(x){
  df <- cellData[cellData$Sample == x,c("nFrags","TSSEnrichment")]
  df$log10nFrags <- log10(df$nFrags)
  df <- df[,c("log10nFrags","TSSEnrichment")]
  identifyCells(df, minTSS=minTSS)
})
names(cellResults) <- samples

# Save models for future reference
saveRDS(cellResults, file = paste0(proj@projectMetadata$outputDirectory, "/cellFiltering.rds"))

# Plot filtering results
for(samp in samples){
  df <- as.data.frame(cellResults[[samp]]$results)
  cell_df <- df[df$classification == "cell",]
  non_cell_df <- df[df$classification != "cell",]
  
  xlims <- c(log10(500), log10(100000))
  ylims <- c(0, 18)
  # QC Fragments by TSS plot w/ filtered cells removed:
  p <- ggPoint(
    x = cell_df[,1], 
    y = cell_df[,2], 
    size = 1.5,
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = xlims,
    ylim = ylims,
    title = sprintf("%s droplets plotted", nrow(cell_df)),
    rastr = TRUE
  )
  # Add grey dots for non-cells
  p <- p + geom_point_rast(data=non_cell_df, aes(x=log10nFrags, y=TSSEnrichment), color="light grey", size=0.5)
  p <- p + geom_hline(yintercept = minTSS, lty = "dashed") + geom_vline(xintercept = log10(1000), lty = "dashed")
  plotPDF(p, name = paste0(samp,"_EM_model_filtered_cells_TSS-vs-Frags.pdf"), ArchRProj = proj, addDOC = FALSE)
}

# Now, filter ATAC project to   only cells
finalCellCalls <- lapply(cellResults, function(x) x$results) %>% do.call(rbind, .)
proj <- addCellColData(proj, data=finalCellCalls$classification, name="cellCall", cells=rownames(finalCellCalls), force=TRUE)
proj <- addCellColData(proj, data=finalCellCalls$cell_uncertainty, name="cellCallUncertainty", cells=rownames(finalCellCalls), force=TRUE)

realCells <- getCellNames(proj)[finalCellCalls$classification == 'cell']
subProj <- subsetArchRProject(proj, cells=realCells, 
                              outputDirectory="filtered_output", dropCells=TRUE, force=TRUE)

# Add sample metadata



# subProj$preservation <- samp.preservation[subProj$Sample2] %>% unlist() %>% as.factor()
# subProj$sex <- samp.sex[subProj$Sample2] %>% unlist() %>% as.factor()
# subProj$age <- samp.age[subProj$Sample2] %>% unlist()

# Now, add tile matrix and gene score matrix to ArchR project
subProj <- addTileMatrix(subProj, force=TRUE)
subProj <- addGeneScoreMatrix(subProj, force=TRUE)

# Add Infered Doublet Scores to ArchR project (~5-10 minutes)
subProj <- addDoubletScores(subProj, dimsToUse=1:20, scaleDims=TRUE, LSIMethod=2)


# Filter doublets:
subProj <- filterDoublets(subProj, filterRatio = 1)

# Save filtered ArchR project
saveArchRProject(subProj)

##########################################################################################
# Reduced Dimensions and Clustering
##########################################################################################

# Add info about whether cell is 'control' or 'disease' (C_SD or AA)
subProj$diseaseStatus <- NA
subProj$diseaseStatus <- ifelse(grepl("C_SD", subProj$Sample), "C_SD", subProj$diseaseStatus)
subProj$diseaseStatus <- ifelse(grepl("C_PB", subProj$Sample), "C_PB", subProj$diseaseStatus)
subProj$diseaseStatus <- ifelse(grepl("AA", subProj$Sample), "AA", subProj$diseaseStatus)



# Reduce Dimensions with Iterative LSI (<5 minutes)
set.seed(1)
subProj <- addIterativeLSI(
  ArchRProj = subProj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  sampleCellsPre = 15000,
  varFeatures = 50000, 
  dimsToUse = 1:25,
  force = TRUE
)

# Identify Clusters from Iterative LSI
subProj <- addClusters(
  input = subProj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.6,
  force = TRUE
)

##########################################################################################
# Visualize Data
##########################################################################################

set.seed(1)
subProj <- addUMAP(
  ArchRProj = subProj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.4, 
  metric = "cosine",
  force = TRUE
)

# Relabel clusters so they are sorted by cluster size

subProj <- addImputeWeights(subProj)


# Save unfiltered ArchR project
saveArchRProject(subProj)

##########################################################################################
# Remove doublet clusters
##########################################################################################

# There are a few clusters that appear to be mostly doublets / poor quality cells
# (Poor TSS enrichment, few marker peaks / GS in downstream analysis, higher cellCallUncertainty)
# and were not filtered by automated doublet removal or by basic QC filtering
# They will be manually filtered here.
proj <- loadArchRProject('filtered_output/', force=TRUE)

# Identify Marker Gene through Pairwise Test vs Bias-Matched Background:
markersGS <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# All valid gene names
geneNames <- rowData(markersGS)$name

# Lists of 'marker genes'
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.00")

# Marker genes we want to highlight for labeling broad clusters:
markerGenes  <- c(
  "KRT1", "KRT10", "KRT5", "KRT14", "KRT15", "DSP", # Keratinocytes
  "THY1", "COL1A1", "COL1A2", "COL3A1", "DCN", "MGP", "COL6A2", # Fibroblasts
  "CD3D", "CD8A", "CD4", "PTPRC", "FOXP3", "IKZF2", "CCL5", # T-cells
  "CD19", "MS4A1", # B-cells
  "CD14", "CD86", "CD74", "CD163", #Monocytes / macrophages
  "VWF", "PECAM1", "SELE", # Endothelial
  "MITF", "TYR", # Melanocyte markers
  "ITGAX", "CD1C", "CD1A", "CLEC1A", "CD207", # Dendritic cells (ITGAX = cd11c)
  "TPM1", "TPM2", "TAGLN" # Muscle
)

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.00", 
  labelMarkers = markerGenes,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  transpose = FALSE
)

hm <- ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(hm, name = "filtered-GeneScores-Marker-Heatmap", width = 6, height = 10, ArchRProj = proj, addDOC = FALSE)

# Remove clusters that have poor quality (enriched for high doublet score cells, poor cell quality, incompatible marker gene scores, etc.)
nonMultipletCells <- getCellNames(proj)[proj$Clusters %ni% c("C7", "C13", "C15", "C18")]

proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = nonMultipletCells,
  outputDirectory = "multiplets_removed_output",
  dropCells=TRUE, force=TRUE
)
saveArchRProject(proj)

# Now, redo clustering and visualization:

# Reduce Dimensions with Iterative LSI
set.seed(1)
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI",
  sampleCellsPre = 20000,
  varFeatures = 50000, 
  dimsToUse = 1:50,
  force = TRUE
)

# Identify Clusters from Iterative LSI
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.7,
  force = TRUE
)

set.seed(1)
proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 60, 
  minDist = 0.6, 
  metric = "cosine",
  force = TRUE
)

# Relabel clusters so they are sorted by cluster size
proj <- relabelClusters(proj)
proj <- addImputeWeights(proj)
# Make various cluster plots:
proj <- visualizeClustering(proj, pointSize=pointSize, sampleCmap=sample_cmap, diseaseCmap=disease_cmap)

# Save filtered ArchR project
saveArchRProject(proj)


##########################################################################################
# Reduced Dimensions and Clustering
##########################################################################################

# Add info about whether cell is 'control' or 'disease' (C_SD or AA)
subProj$diseaseStatus <- NA
subProj$diseaseStatus <- ifelse(grepl("C_SD", subProj$Sample), "C_SD", subProj$diseaseStatus)
subProj$diseaseStatus <- ifelse(grepl("C_PB", subProj$Sample), "C_PB", subProj$diseaseStatus)
subProj$diseaseStatus <- ifelse(grepl("AA", subProj$Sample), "AA", subProj$diseaseStatus)
disease_cmap <- head(cmaps_BOR$stallion,3)
names(disease_cmap) <- c("AA", "C_SD", "C_PB")

# Reduce Dimensions with Iterative LSI (<5 minutes)
set.seed(1)
subProj <- addIterativeLSI(
  ArchRProj = subProj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  sampleCellsPre = 15000,
  varFeatures = 50000, 
  dimsToUse = 1:25,
  force = TRUE
)

# Identify Clusters from Iterative LSI
subProj <- addClusters(
  input = subProj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.6,
  force = TRUE
)

##########################################################################################
# Visualize Data
##########################################################################################

set.seed(1)
subProj <- addUMAP(
  ArchRProj = subProj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.4, 
  metric = "cosine",
  force = TRUE
)

# Relabel clusters so they are sorted by cluster size
subProj <- relabelClusters(subProj)

subProj <- addImputeWeights(subProj)

# Make various cluster plots:
subProj <- visualizeClustering(subProj, pointSize=pointSize, sampleCmap=samp_cmap, diseaseCmap=disease_cmap)

# Save unfiltered ArchR project
saveArchRProject(subProj)

##########################################################################################
# Remove doublet clusters
##########################################################################################

# There are a few clusters that appear to be mostly doublets / poor quality cells
# (Poor TSS enrichment, few marker peaks / GS in downstream analysis, higher cellCallUncertainty)
# and were not filtered by automated doublet removal or by basic QC filtering
# They will be manually filtered here.
proj <- loadArchRProject("/filtered_output/", force=TRUE)

# Identify Marker Gene through Pairwise Test vs Bias-Matched Background:
markersGS <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# All valid gene names
geneNames <- rowData(markersGS)$name

# Lists of 'marker genes'
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.00")

# Marker genes we want to highlight for labeling broad clusters:
markerGenes  <- c(
  "KRT1", "KRT10", "KRT5", "KRT14", "KRT15", "DSP", # Keratinocytes
  "THY1", "COL1A1", "COL1A2", "COL3A1", "DCN", "MGP", "COL6A2", # Fibroblasts
  "CD3D", "CD8A", "CD4", "PTPRC", "FOXP3", "IKZF2", "CCL5", # T-cells
  "CD19", "MS4A1", # B-cells
  "CD14", "CD86", "CD74", "CD163", #Monocytes / macrophages
  "VWF", "PECAM1", "SELE", # Endothelial
  "MITF", "TYR", # Melanocyte markers
  "ITGAX", "CD1C", "CD1A", "CLEC1A", "CD207", # Dendritic cells (ITGAX = cd11c)
  "TPM1", "TPM2", "TAGLN" # Muscle
)

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.00", 
  labelMarkers = markerGenes,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  transpose = FALSE
)

hm <- ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(hm, name = "filtered-GeneScores-Marker-Heatmap", width = 6, height = 10, ArchRProj = proj, addDOC = FALSE)
nonMultipletCells <- getCellNames(proj)[proj$Clusters %ni% c("C6", "C11")]

proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = nonMultipletCells,
  outputDirectory = "multiplets_removed_output",
  dropCells=TRUE, force=TRUE
)
saveArchRProject(proj)
# Reduce Dimensions with Iterative LSI
set.seed(1)
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI",
  sampleCellsPre = 20000,
  varFeatures = 50000, 
  dimsToUse = 1:50,
  force = TRUE
)

# Identify Clusters from Iterative LSI
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.7,
  force = TRUE
)

set.seed(1)
proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 60, 
  minDist = 0.6, 
  metric = "cosine",
  force = TRUE
)

proj <- relabelClusters(proj)
proj <- addImputeWeights(proj)
saveArchRProject(proj)





p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj), 
  plotAs="points"
)
plotPDF(p, name = "Plot-UMAP-marker_gene.pdf", ArchRProj=proj, addDOC=FALSE, width=5, height=5)

labelNew <- c(
  "C1" = "T","C2" = "Keratinocytes","C3" = "B","C4" = "Endothelial","C5" = "Fibroblasts",
  "C6" = "Keratinocytes","C7" = "Fibroblasts","C8" = "Keratinocytes","C9" = "Muscle","C10" = "Keratinocytes",
  "C11" = "T","C12" = "T","C13" = "Keratinocytes","C14" = "Keratinocytes","C15" = "Others",
  "C16" = "Others","C17" = "Monocytes","C18" = "Muscle","C19" = "Keratinocytes","C20" = "B",
  "C21" = "Others","C22" = "B","C23" = "Keratinocytes"
)
proj$cell_type <- mapvalues(proj$Clusters, names(labelNew), labelNew)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "cell_type")
plotPDF(p1, name = "Plot-UMAP-cell_type.pdf", ArchRProj=proj, addDOC=FALSE, width=5, height=5)

final_cell <- getCellNames(proj)[proj$cell_type %ni% c("Others")]
proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = final_cell,
  outputDirectory = "multiplets_removed_output",
  dropCells=TRUE, force=TRUE
)
saveArchRProject(proj)   
proj <- addGroupCoverages(
  ArchRProj=proj, 
  groupBy="cell_type", 
  minCells = 50, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40)
  force=TRUE
)
proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "cell_type", 
  peaksPerCell = 500, # The upper limit of the number of peaks that can be identified per cell-grouping in groupBy. (Default = 500)
  pathToMacs2 = '/home/kyh/miniconda3/envs/RNA-seq/envs/deeptools/bin/macs2',
  force = TRUE
)
proj <- addPeakMatrix(ArchRProj = proj, force = TRUE)

saveArchRProject(proj)

sample28_peak_matrix <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE
)
saveRDS(sample28_peak_matrix,file = 'Sample28_peak_matrix.Rds')

subset_peaks <- function(x){
  if (length(str_extract_all(x,'-')[[1]])== 2) {
    if (str_starts(x,'chr')) {
      return(TRUE)
    }else{
      return(FALSE)
    }
    return(TRUE)
  }else{
    return(FALSE)
  }
}

sparse_mtx <- sample28_peak_matrix@assays@data$PeakMatrix
rownames(sparse_mtx) <- paste(as.data.frame(sample28_peak_matrix@rowRanges)[[1]],as.data.frame(sample28_peak_matrix@rowRanges)[[2]],as.data.frame(sample28_peak_matrix@rowRanges)[[3]],sep = '-')
sparse_mtx <- sparse_mtx[which(purrr::map_vec(rownames(sparse_mtx),subset_peaks)),]
library(Seurat)
library(Signac)
chrom_assay <- CreateChromatinAssay(
  counts = sparse_mtx,
  sep = c("-", "-"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)
sample28_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = as.data.frame(sample28_peak_matrix@colData)
)
saveRDS(sample28_ATAC,file = 'sample28_scATAC-seq_29k_processed.Rds')
