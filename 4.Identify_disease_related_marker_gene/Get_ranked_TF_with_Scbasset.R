library(Signac)
library(Seurat)
library(tidyverse)
library(Matrix)
library(glue)
source('../2.Data preprocessing//preprocess_functions.R')


get_celltype_TF_act <- function(df){
    cell_type_tf_act <- map_dfr(1:nrow(df),function(i){
        cell_type <-  df[i,1] %>% varhandle::unfactor() %>% as.character() 
        cell_type_top10 <- df[i,-1] |>  t() |> as.data.frame() %>% 
            rownames_to_column('gene') %>% slice_max(order_by = V1,n = 10) %>% 
            mutate(cell_types = cell_type) %>% 
            dplyr::select(cell_types,gene)
        return(cell_type_top10)
    },.progress = T)
    return(cell_type_tf_act)
}


##
sample_rds <- list.files('../../data/Rds/',full.names = T)
names(sample_rds) <- str_extract(sample_rds,'sample\\d+|Sample\\d+')

tf_activity <- list.files('../../data/TF_activity/rawdata/',full.names = T)
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
    
    out_file <- paste0('../../data/TF_activity/Temp/',names(sample_rds)[i],'cell_type_activity.Rds')
    saveRDS(TF_act_metadata,out_file)
    gc()
    cat(names(sample_rds)[i],'is',class(rds),'Object finished\n')
}


sample12 <- readRDS('../../data/Rds/Sample12_ATAC-seq.Rds')
sample12_TF_act_metadata <- sample12@meta.data %>% rownames_to_column('cell_barcode')
sample12_TF_act <- data.table::fread('../../data/TF_activity/rawdata/sample12_scATAC-seq_104k_processed_tf_activity.txt.gz')
sample12_TF_act_metadata <- left_join(sample12_TF_act_metadata,sample12_TF_act,by = c('cell_barcode' = 'V1')) %>% 
    filter(!is.na(cell_type)) %>%
    group_by(cell_type) %>% 
    summarise(across(where(is.numeric),sum)) %>% 
    dplyr::select(c('cell_type',colnames(sample12_TF_act)[-1]))


sample12_TF_act_metadata$cell_type <- str_replace_all(sample12_TF_act_metadata$cell_type,'/','')

##
sample5_rds <- list.files('../../data/scATAC-seq/Sample5/',pattern = 'RDS',full.names = T)
sample5_TF_act <- data.table::fread('../../data/TF_activity/rawdata/sample5_scATAC-seq_720k_processed_tf_activity.txt')
sample5_metadata <- map(sample5_rds,function(x){
    sce <- readRDS(x)
    metadata <- sce@meta.data %>% rownames_to_column('cell_barcode')
    return(metadata)
},.progress = T)
sample5_metadata_df <- bind_rows(sample5_metadata)
sample5_TF_act_metadata <- left_join(sample5_metadata_df,sample5_TF_act,by = c('cell_barcode' = 'V1')) %>% 
    filter(!is.na(cell_type)) %>%
    group_by(cell_type) %>% 
    summarise(across(where(is.numeric),sum)) %>% 
    dplyr::select(c('cell_type',colnames(sample5_TF_act)[-1]))
saveRDS(sample5_TF_act_metadata,'../../data/TF_activity/Temp/sample5cell_type_activity.Rds')

###

sample11_rds <- readRDS(sample_rds[3])
sample11_TF_act_metadata <- sample11_rds@meta.data %>% rownames_to_column('cell_barcode')
sample11_TF_act <- data.table::fread('../../data/TF_activity/rawdata/sample11_scATAC-seq_17k_processed_tf_activity.txt',sep = '\t')


sample11_TF_act_metadata <-
    left_join(sample11_TF_act_metadata,
              sample11_TF_act,
              by = c('cell_barcode' = 'V1')) %>%
    filter(!is.na(cell_type))  %>%
    group_by(cell_type) %>%
    summarise(across(where(is.numeric), sum)) %>%
    dplyr::select(c('cell_type', colnames(sample11_TF_act)[-1]))
saveRDS(sample11_TF_act_metadata,'../../data/TF_activity/Temp/sample11cell_type_activity.Rds')

Sample11_tf_activity <- get_celltype_TF_act(sample11_TF_act_metadata)
data.table::fwrite(Sample11_tf_activity,
                   file =  '../../data/TF_activity/sample11_tf_activity_top10.txt',
                   sep = '\t')






###
tf_activity_rds <- list.files('../../data/TF_activity/Temp/',pattern = 'Rds',full.names = T)
tf_activity_rds_list <- map(tf_activity_rds,function(x){
    tf_activity <- readRDS(x)
    tf_activity_res <- get_celltype_TF_act(tf_activity)
    return(tf_activity_res)
}) %>% setNames(str_extract(tf_activity_rds,'sample\\d+|Sample\\d+'))
saveRDS(tf_activity_rds_list,'../../data/TF_activity/tf_activity_rds_list.Rds')


tf_activity_rds_list <- readRDS('../../data/TF_activity/tf_activity_rds_list.Rds')
for (i in 1:length(tf_activity_rds_list)) {
    filename <- paste0('../../data/TF_activity/',names(tf_activity_rds_list)[i],'_tf_activity_top10.txt')
    data.table::fwrite(tf_activity_rds_list[[i]],file = filename,sep = '\t')
}

Sample12_tf_activity <- get_celltype_TF_act(sample12_TF_act_metadata)

Sample12_tf_activity <- as.data.frame(Sample12_tf_activity)

Sample12_tf_activity2 <-
    tibble(cell_types = map_chr(1:nrow(Sample12_tf_activity), function(i) {
        Sample12_tf_activity$cell_types[i, 1] |> as.character()
    }),
    gene = Sample12_tf_activity$gene)


data.table::fwrite(Sample12_tf_activity2,
            file =  paste0('../../data/TF_activity/',names(tf_activity_rds_list)[2],'_tf_activity_top10.txt'),
            sep = '\t')

Sample5_tf_activity <- get_celltype_TF_act(sample5_TF_act_metadata)
data.table::fwrite(Sample12_tf_activity2,
                   file =  '../../data/TF_activity/sample5_tf_activity_top10.txt',
                   sep = '\t')


