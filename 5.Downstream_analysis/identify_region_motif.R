peak_anno <- function(merged_footprints, tssRegion = c(-3000, 3000)) {
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  library(org.Hs.eg.db)
  merged_footprints[,2] = as.numeric(merged_footprints[,2])
  merged_footprints[,3] = as.numeric(merged_footprints[,3])
  annodb <- 'org.Hs.eg.db'
  
  reference_GRange <- GenomicRanges::GRanges(seqnames = merged_footprints[,1],
                                             IRanges::IRanges(start = as.numeric(merged_footprints[,2]),
                                                              end = as.numeric(merged_footprints[,3])))
  peakAnno <- ChIPseeker::annotatePeak(reference_GRange,
                                       tssRegion = tssRegion,
                                       TxDb = txdb, annoDb = annodb
  )
  region <- peakAnno@anno@elementMetadata$annotation
  gene <- peakAnno@anno@elementMetadata$ENSEMBL
  symbol <- peakAnno@anno@elementMetadata$SYMBOL
  start1 <- peakAnno@anno@ranges@start
  dis <- peakAnno@anno@elementMetadata$distanceToTSS
  merged_footprints2 <- merged_footprints[merged_footprints$V2 %in% start1, ]
  exon1 <- grep('exon',region)
  Intron1 <- grep('Intron',region)
  Intergenic1 <- grep('Intergenic',region)
  Downstream1 <- grep('Downstream',region)
  Promoter1 <- grep('Promoter',region)
  UTR3 <- grep("3' UTR",region)
  UTR5 <- grep("5' UTR",region)
  region2 <- rep(NA,length(region))
  region2[exon1]='Exon'
  region2[Intron1]='Intron'
  region2[Downstream1]='Downstream'
  region2[Promoter1]='Promoter'
  region2[UTR3]="3' UTR"
  region2[UTR5]="5' UTR"
  region2[Intergenic1]='Intergenic'
  table(region2)
  peak_region1 <- paste(as.character(peakAnno@anno@seqnames),
                        as.character(peakAnno@anno@ranges),sep = ':')
  peak_region2 <- paste0(merged_footprints[,1],':'
                         ,merged_footprints[,2],'-',merged_footprints[,3])
  merged_footprints2 <- merged_footprints[peak_region2%in%peak_region1,]
  peak_gr=GenomicRanges::GRanges(seqnames = merged_footprints2[,1],
                                 IRanges::IRanges(start=merged_footprints2[,2],
                                                  end=merged_footprints2[,3]))
  peak_gr$gene = gene
  peak_gr$symbol = symbol
  peak_gr$region = region2
  peak_gr$distanceToTSS = dis
  
  ### filter peaks based on annotation
  peak_gr = peak_gr[peak_gr$region=='Promoter']
  return(peak_gr)
}

Str_to_GR <- function(x){
  sp = strsplit(x,split='-')
  chr = sapply(sp,function(x) x[[1]])
  s = as.numeric(sapply(sp,function(x) x[[2]]))
  e = as.numeric(sapply(sp,function(x) x[[3]]))
  GR_out = GRanges(chr,IRanges(s,e))
  return(GR_out)
}

Seqnames <- function(x){
  out= as.character(seqnames(x))
  return(out)
}

Start <- function(x){
  out= as.numeric(start(x))
  return(out)
}

End <- function(x){
  out= as.numeric(end(x))
  return(out)
}

Must_to_GR <- function(x){
  library('pbapply')
  chr_all = pblapply(x,Seqnames)
  start_all = pblapply(x,Start)
  end_all = pblapply(x,End)
  len_all = pblapply(x,function(x) length(x))
  #####
  chr_all = as.character(unlist(chr_all))
  start_all = as.numeric(unlist(start_all))
  end_all = as.numeric(unlist(end_all))
  ##### ##### 
  names_all = rep(names(x),len_all)
  GR_out = GRanges(chr_all,IRanges(start_all,end_all),motifs=names_all)
  return(GR_out)
}

motifs_select <- function(motif,gene){
  index <- c()
  if (stringr::str_sub(gene[1],1,3)=='ENS') {
    col_idx = 5
  }else{col_idx = 4}
  for (i in 1:nrow(motif)) {
    judge <- c()
    gene1 <- strsplit(motif[i,col_idx],';')[[1]]
    for (j in gene1) {
      if (j %in% gene) {
        judge <- c(judge,'YSE')
      }
    }
    if ('YSE' %in% judge) {
      index <- c(index,i)
    }
  }
  motif1 <- motif[index,]
  return(motif1)
}

identify_region_tfs <- function(All_peaks_GR,gene.use,PWM,motifdb,pvalue.cutoff = 5e-05){
  library(motifmatchr)
  library(BSgenome.Hsapiens.UCSC.hg38)
  motif_use = motifs_select(motifdb,gene.use)
  PWM = PWM[motif_use$Accession]
  matched_motif <- matchMotifs(PWM,
              All_peaks_GR,genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
              out='positions',p.cutoff = pvalue.cutoff)
  matched_motif <- Must_to_GR(matched_motif)
  return(matched_motif)
}

overlap_peak_motif <- function(peak,motif,motifdb){
  overlaped = findOverlaps(peak,motif)
  peak_motif = cbind(as.data.frame(peak[overlaped@from]),as.data.frame(motif[overlaped@to]))
  peak_motif$TF = motifdb[match(peak_motif$motifs,motifdb$Accession),4]
  return(peak_motif)
}

make_tf_target <- function(atac_out){
  tf = atac_out$TF
  tf = paste0(tf,'#',atac_out$symbol)
  tf_target = unlist(map(tf,~paste_gene(.x)))
  return(tf_target)
}

paste_gene <- function(gene){
  tf = strsplit(gene,'#')[[1]][1]
  target = strsplit(gene,'#')[[1]][2]
  tf = unlist(strsplit(tf,';'))
  return(paste0(tf,'-',target))
}
