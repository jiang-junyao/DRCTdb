library(tidyverse)
library(CSCORE)
library(Seurat)
# g <- graph_from_data_frame(test2, directed=FALSE)
# V(g)$size <- 5
# pdf('test.pdf',width = 10,height = 10)
# plot(g, vertex.size=V(g)$size, edge.width=E(g)$width, vertex.label.cex=0.3)
# dev.off()

get_cor <- function(df,sce){
    grn <- df %>% select(TF,symbol) %>%separate_rows(TF,sep = ';') 
    if ('nCount_RNA' %in% sce@meta.data) {
        mean_exp <-  rowMeans(sce@assays$RNA@counts/sce$nCount_RNA)
    }else{
        mean_exp <-  rowMeans(sce@assays$RNA@counts/sce$nCounts_RNA)
    }
   
    mean_exp <- names(mean_exp[which(mean_exp > 0.1) ])
    genes_selected  <-  c(unique(grn$TF),unique(grn$symbol)) %>% 
        intersect(rownames(sce)) %>%  intersect(mean_exp) 
    beta_CSCORE_result <- CSCORE(sce, genes = genes_selected)
    
    
    # Obtain CS-CORE co-expression estimates
    CSCORE_coexp <- beta_CSCORE_result$est
    # Obtain BH-adjusted p values
    CSCORE_p <- beta_CSCORE_result$p_value
    p_matrix_BH = matrix(0, length(genes_selected), length(genes_selected))
    p_matrix_BH[upper.tri(p_matrix_BH)] = p.adjust(CSCORE_p[upper.tri(CSCORE_p)], method = "BH")
    p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)
    
    # Set co-expression entires with BH-adjusted p-values greater than 0.05 to 0
    CSCORE_coexp[p_matrix_BH > 0.05] <- 0
    grn <- grn %>% filter(TF %in% genes_selected, symbol %in% genes_selected) 
    grn$cor <- map_dbl(1:nrow(grn),function(i){
        CSCORE_coexp[grn$TF[i],grn$symbol[i]]
    },.progress = T)
    
    refined_grn <- grn %>% filter(cor != 0) %>% filter(abs(cor) > 0.3,TF != symbol ) %>% distinct_all()
    gc()
    return(refined_grn)
    
}
batch_run_cor <- function(file,sce){
    if (!('Seurat' %in% class(sce))) {
        stop('Must input Seurat object')
    }
    if (!('cell_type' %in% colnames(sce@meta.data))) {
        stop('Must have cell_type in metadata')
    }
    
    filename <-  tools::file_path_sans_ext(basename(file))
    sample <- str_extract(file,'sample\\d+')
    cell_type <- str_extract(filename,'.*(?=#)')
    subset_sce <- subset(sce,cell_type == cell_type)
    
    tf_target <- readRDS(file)
    tf_target_with_cor <- get_cor(tf_target,subset_sce)
    out_name <- paste0('../../data/tf_target/',sample,'/',filename,'with_cor.RDS')
    saveRDS(tf_target_with_cor,out_name)
    cat(out_name,'finished\n')
    return(tf_target_with_cor)
    
}

tf_target <- readRDS('../../data/downstream_result/sample23/pval0001_tftarget.rds')
sample23_islet_scRNA <- readRDS('../../data/Rds/RNA/sample23_islet_scRNA_95k_processed.Rds')

beta <- subset(sample23_islet_scRNA,cell_type == 'Beta') 
beta_GRN <- get_cor(tf_target$Beta_Type_2_diabetes,beta)
Delta <- subset(sample23_islet_scRNA,cell_type == 'Delta') 
Delta_GRN <- get_cor(tf_target$Delta_Type_2_diabetes,Delta)
tf_target$Beta_Type_2_diabetes <- beta_GRN
tf_target$Delta_Type_2_diabetes <-Delta_GRN
saveRDS(tf_target,file = '../../data/GRN/pval0001_tftarget_with_cor.rds')


total_rds <- list.files('../../data/tf_target/',recursive = T,full.names = T)

#sample1----
sample1 <- readRDS('../../data/Rds/RNA/sample1_heart_scRNA_35k_processed.Rds')
sample1$cell_type <- sample1$celltype
sample1_tf_core <- str_subset(total_rds,'sample1/')

sample1_cor <- map(
    sample1_tf_core,batch_run_cor,sce = sample1,.progress = T
)
saveRDS(sample1_cor,'../../data/GRN/sample1_GRN.Rds')
#sample2----
sample3 <- readRDS('../../data/Rds/RNA/sample3_atlas_scRNA_600k_processed.Rds')
sample3_tf_core <- str_subset(total_rds,'sample3/')
sample3_cor <- map(
    sample3_tf_core,batch_run_cor,sce = sample3,.progress = T
)
saveRDS(sample3_cor,'../../data/GRN/sample3_GRN.Rds')
