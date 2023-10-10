library(tidyverse)

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


ldsc_res <- list.files('../../data/LDSC_results/',pattern = 'pvalues.tsv',recursive = T)
    
ldsc_res_list <- map(paste0('../../data/LDSC_results/',ldsc_res),function(x){
    df <- data.table::fread(x) %>% as_tibble() %>% 
        pivot_longer(cols = -V1,names_to = 'cell_type',values_to = 'pvalue') %>% 
        mutate(neglog10p = -log10(pvalue)) %>% 
        arrange(desc(neglog10p)) %>%
        filter(abs(neglog10p) >= 3)  %>% 
        select(1:3) 
    colnames(df) <- c('Disease', 'Cell type','p value')
    return(df)
}) %>% setNames(str_extract(ldsc_res,'\\w+'))


writexl::write_xlsx(ldsc_res_list,'../../data/LDSC_results/DRCTS_res.xlsx') 

DRCTS <- map(ldsc_res_list,function(x){
    unique(x[['Disease']])
})


DRCT_stat <- tibble(
    dataset = names(DRCTS),
    num = map_int(DRCTS,length),
    order = as.numeric(str_extract(dataset,'\\d+'))
) %>% filter(dataset != '') %>% arrange(order)
clustcol<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")


DRCT_stat$dataset <- factor(DRCT_stat$dataset,levels = DRCT_stat$dataset)
p1 <- ggplot(DRCT_stat,aes(dataset,num,fill = dataset)) +
    labs(x = 'dataset', y = 'Enriched diease numbers') +
    geom_bar(stat = "identity",position = position_dodge(width = 1)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(labels = c(1:15)) +
    scale_fill_manual(values = clustcol[1:nrow(DRCT_stat)]) +
    theme_bw() +
    guides(fill = F) 
p1
cairo_pdf('../../data/enriched_disease2.pdf',width = 4,height = 2.5)
p1
dev.off()

ggplot(test_res,aes(gene_name,Fold,fill= sample)) + 
    geom_bar(stat = "identity",position = position_dodge(width = 1)) + theme_classic() +
    geom_errorbar(aes(ymin = Fold - SD, ymax = Fold + SD),
                  position=position_dodge(width=1), 
                  width=0.3,size=0.3,colour="black") +
    labs(x = '',y = 'Relative expression levels')+
    guides()+
    ggprism::theme_prism()+
    scale_fill_manual(values = c("#69b3a2","#836FFF"))+
    theme(axis.title.x = element_blank(),
          aspect.ratio = 1.3,
          axis.text.x = element_text(angle = 45))+
    ggpubr::geom_signif(y_position = get_maxfold(test_res), 
                        xmin = rep(1:length(get_ano(test_res))) -0.2, 
                        xmax = rep(1:length(get_ano(test_res))) + 0.2,
                        annotation = get_ano(test_res),
                        tip_length = 0)


DRCTS_res <- tibble(
    dataset = names(DRCTS),
    disease = map_chr(DRCTS,function(x){
        if(length(x) > 3){
            x <- x[1:3]
            x <- paste0(x,collapse = ';') %>% paste0('...')
        }else{
            x <- paste0(x,collapse = ';')
        }
        
        return(x)
    })
) %>% arrange(as.numeric(str_extract(dataset,'\\d+')))

data.table::fwrite(DRCTS_res,file = '../../data/LDSC_results/DRCTS_res.xls',sep = '\t')






