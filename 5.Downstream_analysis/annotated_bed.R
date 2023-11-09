library(ChIPseeker)
library(tidyverse)
library(GenomicRanges)
peak_anno <- function(reference_GRange, tssRegion = c(-2000, 500),filter = FALSE) {
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    library(org.Hs.eg.db)
    annodb <- 'org.Hs.eg.db'
    
    peakAnno <- ChIPseeker::annotatePeak(reference_GRange,
                                         tssRegion = tssRegion,
                                         TxDb = txdb, annoDb = annodb,
                                         verbose = F
    )
    region <- peakAnno@anno@elementMetadata$annotation
    gene <- peakAnno@anno@elementMetadata$ENSEMBL
    symbol <- peakAnno@anno@elementMetadata$SYMBOL
    start1 <- peakAnno@anno@ranges@start
    dis <- peakAnno@anno@elementMetadata$distanceToTSS
    
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
    peak_gr=reference_GRange
    peak_gr$gene = gene
    peak_gr$symbol = symbol
    peak_gr$region = region2
    peak_gr$distanceToTSS = dis
    ### filter peaks based on annotation
    if (filter) {
        peak_gr = peak_gr[peak_gr$region=='Promoter']  
    }
    return(peak_gr)
}

all_DERs <- list.files('../../data/DERs/',pattern = 'txt$',full.names = T)
map(all_DERs,function(file){
    df <-
        data.table::fread(file) %>% separate(
            col = names,
            into = c('seqnames', 'start', 'end'),
            sep = '-'
        ) |>
        dplyr::select(seqnames, start, end, group, logfoldchanges, pvals_adj) |> 
        makeGRangesFromDataFrame(keep.extra.columns = T) |> 
        peak_anno() |>
        as.data.frame()  
    data.table::fwrite(df,file = file,sep = '\t')
    cat(file,'Finished\n')
},.progress = T)

for (i in list.files('../../data/DERs/')){
    raw_path <- paste0('../../data/DERs/',i)
    new_path <- paste0('../../data/downstream_result/',str_extract(i,'sample\\d+'),'/',i)
    file.copy(raw_path,new_path,overwrite = T)
}





snp_file <- list.files('../../data/downstream_result/',pattern = '.txt',full.names = T,recursive = T) %>% str_subset('snp')

map(snp_file,function(file){
    df <-
        data.table::fread(file) 
    colnames(df)[7] <- 'SNP_location'
    gr <- df %>% dplyr::select(c(1:3, 11, 7)) %>% dplyr::mutate(seqnames = paste0('chr',seqnames)) %>%
        makeGRangesFromDataFrame(keep.extra.columns = T) 
    df1 <-  peak_anno(gr) %>% as.data.frame()  
    data.table::fwrite(df1,file = file,sep = '\t')
},.progress = T)

snp_file <- list.files('../../data/downstream_result/sample28/',pattern = '.txt',full.names = T,recursive = T) %>% str_subset('snp')
map(snp_file,function(file){
    df <-
        data.table::fread(file) 
    colnames(df)[7] <- 'SNP_location'
    gr <- df %>% dplyr::select(c(1:3, 11, 7)) %>% dplyr::mutate(seqnames = paste0('chr',seqnames)) %>%
        makeGRangesFromDataFrame(keep.extra.columns = T) 
    df1 <-  peak_anno(gr) %>% as.data.frame()  
    data.table::fwrite(df1,file = file,sep = '\t')
},.progress = T)


##reduce the file size of DERs
raw_path = list.files('../../data/DERs/',full.names = T)
DER_list <- map(raw_path,data.table::fread)


DER_list_small <- map(DER_list,function(df){
    if (nrow(df) >1000) {
        df <- df %>% group_by(seqnames) %>% arrange(desc(abs(logfoldchanges))) %>%
            slice_head(n = 3000) 
    }
    df$logfoldchanges <-  round(df$logfoldchanges,2)
    df$pvals_adj <- round(df$pvals_adj,3)
    return(df)
})
for (i in 1:length(DER_list_small)) {
    data.table::fwrite(DER_list_small[[i]],file = raw_path[i],sep = '\t')
}
