library(Signac)
library(Seurat)
library(tidyverse)
library(Matrix)
library(glue)
source('preprocess_functions.R')

sample_rds <- list.files('../../data/Rds/',full.names = T)[-3]
names(sample_rds) <- str_extract(sample_rds,'sample\\d+|Sample\\d+')

tf_activity <- list.files('../../data/TF_activity/',full.names = T)[-12]
names(tf_activity) <- str_extract(tf_activity,'sample\\d+|Sample\\d+')




for (i in 1:length(sample_rds)) {
    rds <- readRDS(sample_rds[i])
    tf_activity_file <- str_subset(tf_activity,fixed(paste0(names(sample_rds)[i],'_'),ignore_case = T) )

    TF_act_metadata <- rds@meta.data %>% rownames_to_column('cell_barcode')
    TF_act <- data.table::fread(tf_activity_file)
    TF_act_metadata <-
        left_join(TF_act_metadata,
                  TF_act,
                  by = c('cell_barcode' = 'V1')) %>%
        filter(!is.na(cell_type)) %>%
        group_by(cell_type) %>%
        summarise(across(where(is.numeric), sum)) %>%
        dplyr::select(c('cell_type', colnames(TF_act)[-1]))
    
    out_file <- paste0('../../data/TF_activity/',names(sample_rds)[i],'cell_type_activity.Rds')
    saveRDS(TF_act_metadata,out_file)
    cat(names(sample_rds)[i],'is',class(rds),'Object finished\n')
}


sample1_ATAC <- readRDS('../../data/scATAC-seq/sample1/Rds/sample1_scATAC-seq_80k_processed.Rds')
sample1_TF_act_metadata <- sample1_ATAC@meta.data %>% rownames_to_column('cell_barcode')
sample1_TF_act <- data.table::fread('../../data/TF_activity/sample1_scATAC-seq_80k_processed_tf_activity.txt.gz')
sample1_TF_act_metadata <- left_join(sample1_TF_act_metadata,sample1_TF_act,by = c('cell_barcode' = 'V1')) %>% 
    filter(!is.na(cell_type)) %>%
    group_by(cell_type) %>% 
    summarise(across(where(is.numeric),sum)) %>% 
    dplyr::select(c('cell_type',colnames(sample1_TF_act)[-1]))
saveRDS(sample1_TF_act_metadata,'../../data/TF_activity/sample1_cell_type_activity.Rds')





get_celltype_TF_act <- function(df){
    cell_type_tf_act <- map_dfr(1:nrow(df),function(i){
        cell_type <-  df[i,1] %>% as.character()
        cell_type_top10 <- df[i,-1] |>  t() |> as.data.frame() %>% 
            rownames_to_column('gene') %>% slice_max(order_by = V1,n = 10) %>% 
            mutate(cell_types = cell_type) %>% 
            dplyr::select(cell_types,gene)
        return(cell_type_top10)
    },.progress = T)
    return(cell_type_tf_act)
}

tf_activity_rds <- list.files('../../data/TF_activity/',pattern = 'Rds',full.names = T)
tf_activity_rds_list <- map(tf_activity_rds,function(x){
    tf_activity <- readRDS(x)
    tf_activity_res <- get_celltype_TF_act(tf_activity)
    return(tf_activity_res)
})








