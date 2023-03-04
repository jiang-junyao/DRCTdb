cal_ct_score <- function(overlap_DBRs,DBRs_gr){
  
  ### process
  overlap_peak <- as.data.frame(overlap_DBRs)
  overlap_peak$peak_name <- paste0(overlap_peak[,1],
                                   ':',overlap_peak[,2],'-',overlap_peak[,3])
  DBRs <- as.data.frame(DBRs_gr)
  DBRs$peak_name <- paste0(DBRs[,1],
                          ':',DBRs[,2],'-',DBRs[,3])
  cluster <- levels(as.factor(DBRs$cluster))
  
  ### Hypergeometric test
  prob <- c()
  for (i in cluster) {
    
    DBRs.use = DBRs[DBRs$cluster %in% i,]
    overlap_peak.use = overlap_peak[overlap_peak$cluster %in% i,]
  
    k <- nrow(overlap_peak.use)
    K <- nrow(overlap_peak)
    N <- nrow(DBRs)
    n <- nrow(DBRs.use)
    
    if (n==0 | is.null(n)) {
      p1 = 1
    }else{
     p1 = phyper(n,k,K,N,lower.tail = FALSE) 
    }
    prob <- c(prob,p1)
  }
  pfdr <- p.adjust(prob,method = 'fdr')
  ct_score <- data.frame(cluster,prob,pfdr)
  colnames(ct_score) <- c('cell_type','pval','fdr')
  
  return(ct_score)
}