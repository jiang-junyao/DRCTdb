
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

identify_region_tfs <- function(All_peaks_GR,pvalue.cutoff = 5e-05){
  library(motifmatchr)
  library(BSgenome.Hsapiens.UCSC.hg38)
  matched_motif <- matchMotifs(PWM,
              All_peaks_GR,genome = BSgenome.Hsapiens.UCSC.hg38.masked,
              out='positions',p.cutoff = pvalue.cutoff)
  matched_motif <- Must_to_GR(Total_footprint_Motif)
  return(matched_motif)
}



identify_regulation <- function(){
  
}