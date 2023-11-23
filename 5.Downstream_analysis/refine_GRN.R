library(tidyverse)
library(CSCORE)
library(Seurat)
tf_target <- readRDS('../../data/downstream_result/sample23/pval0001_tftarget.rds')

test <- tf_target$Beta_Type_2_diabetes %>% select(TF,symbol) %>%
    separate_rows(TF,sep = ';') 

sample23_islet_scRNA <- readRDS('../../data/Rds/sample23_islet_scRNA_95k_processed.Rds')


beta <-  sample23_islet_scRNA[,sample23_islet_scRNA$cell_type %in% 'Beta']
mean_exp <-  rowMeans(beta@assays$RNA@counts/beta$nCount_RNA)
mean_exp <- names(mean_exp[which(mean_exp > 0) ])

genes_selected  <-  c(unique(test$TF),unique(test$symbol)) %>% 
    intersect(rownames(beta)) %>%  intersect(mean_exp) 

beta_CSCORE_result <- CSCORE(beta, genes = genes_selected)



# Obtain CS-CORE co-expression estimates
CSCORE_coexp <- beta_CSCORE_result$est
# Obtain BH-adjusted p values
CSCORE_p <- beta_CSCORE_result$p_value
p_matrix_BH = matrix(0, length(genes_selected), length(genes_selected))
p_matrix_BH[upper.tri(p_matrix_BH)] = p.adjust(CSCORE_p[upper.tri(CSCORE_p)], method = "BH")
p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)

# Set co-expression entires with BH-adjusted p-values greater than 0.05 to 0
CSCORE_coexp[p_matrix_BH > 0.05] <- 0
test <- test %>% filter(TF %in% genes_selected, symbol %in% genes_selected) 
test$cor <- map_dbl(1:nrow(test),function(i){
    CSCORE_coexp[test$TF[i],test$symbol[i]]
},.progress = T)

test2 <- test %>% filter(cor != 0) %>% filter(abs(cor) > 0.3,TF != symbol ) %>% distinct_all()





g <- graph_from_data_frame(test2, directed=FALSE)
V(g)$size <- 5
pdf('test.pdf',width = 10,height = 10)
plot(g, vertex.size=V(g)$size, edge.width=E(g)$width, vertex.label.cex=0.3)
dev.off()



get_cor <- function(df,sce){
    grn <- df %>% select(TF,symbol) %>%separate_rows(TF,sep = ';') 
    mean_exp <-  rowMeans(sce@assays$RNA@counts/sce$nCount_RNA)
    mean_exp <- names(mean_exp[which(mean_exp > 0) ])
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
    return(refined_grn)
    
}


tf_target <- readRDS('../../data/downstream_result/sample23/pval0001_tftarget.rds')
sample23_islet_scRNA <- readRDS('../../data/Rds/sample23_islet_scRNA_95k_processed.Rds')


beta <- subset(sample23_islet_scRNA,cell_type == 'Beta') 
beta_GRN <- get_cor(tf_target$Beta_Type_2_diabetes,beta)

Delta <- subset(sample23_islet_scRNA,cell_type == 'Delta') 
Delta_GRN <- get_cor(tf_target$Delta_Type_2_diabetes,Delta)

tf_target$Beta_Type_2_diabetes <- beta_GRN
tf_target$Delta_Type_2_diabetes <-Delta_GRN
saveRDS(tf_target,file = '../../data/GRN/pval0001_tftarget_with_cor.rds')
