library(tidyverse)
library(data.table)

drct_file <- list.files('./downstream_result/',recursive = T,pattern = 'related_disease.xls',full.names = T)

map(drct_file,function(files){
    drct <- fread(files)
    row_index <- which(str_detect(drct$Disease,'Myocardial fractal dimension'))
    if (length(row_index) > 0) {
        drct$Disease[row_index] <- stringr::str_extract(drct$Disease[row_index],'[\\s\\w]+')
        fwrite(drct,file = files,sep = '\t')
    }
},.progress = T)


all_drct <- map_dfr(drct_file,function(files){
    drct <- fread(files)
    row_index <- which(str_detect(drct$`Cell type`,'FB2$'))
    if (length(row_index) > 0) {
        drct <- drct[-row_index,]
        fwrite(drct,file = files,sep = '\t')
        return(drct)
    }else{
        fwrite(drct,file = files,sep = '\t')
        return(drct)
    }
})




all_Myocardial_fractal <- list.files('downstream_result/',recursive = T,pattern = 'Myocardial_fractal_dimension',full.names = T)
for (file in all_Myocardial_fractal) {
    dir_name <- dirname(file)
    file_suffix <- tools::file_ext(file)
    base_name <- basename(file)
    new_base_name <- str_extract(base_name,'\\w+')
    new_name <- paste0(dir_name,'/',new_base_name,'.',file_suffix)
    file.rename(from = file,to = new_name)
}

#
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

total_atac_file <- list.files('./downstream_result',pattern = 'txt',recursive = T,full.names = T) %>%
    str_subset('atac_snp')

atac_table <- data.table::fread(total_atac_file[1])
atac_table_enrich_res <- enrichGO(atac_table$symbol,
                                  keyType = 'SYMBOL',
                                  OrgDb = org.Hs.eg.db,
                                  ont = 'ALL',
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05)
p1 <- dotplot(atac_table_enrich_res,label_format = 40)
saveplot(p1,'Fig4B',width = 4,height = 2)

map(total_atac_file,function(file_path){
    filename <- tools::file_path_sans_ext(file_path)
    atac_table <- data.table::fread(paste0(filename,'.txt'))
    atac_table_enrich_res <- enrichGO(atac_table$symbol,
                                      keyType = 'SYMBOL',
                                      OrgDb = org.Hs.eg.db,
                                      ont = 'ALL',
                                      pvalueCutoff = 0.05,
                                      qvalueCutoff = 0.05)
    p1 <- dotplot(atac_table_enrich_res)
    ggsave(filename = paste0(filename,'.png'),p1,width = 200,height = 100,units = 'mm',dpi = 300) 
},.progress = T)
total_rna_file <- list.files('./downstream_result',pattern = 'txt',recursive = T,full.names = T) %>%
    str_subset('rna_snp')

rna_table <- data.table::fread(total_rna_file[1])
rna_table_enrich_res <- enrichGO(rna_table$symbol,
                                  keyType = 'SYMBOL',
                                  OrgDb = org.Hs.eg.db,
                                  ont = 'ALL',
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05)
p2 <- dotplot(rna_table_enrich_res,label_format = 50)
saveplot(p2,'Fig4C',width = 4,height = 2)


rna_table <- data.table::fread('downstream_result/sample23/rna_snp/Beta_Type_2_diabetes.txt')
rna_table_enrich_res <- enrichGO(rna_table$symbol,
                                 keyType = 'SYMBOL',
                                 OrgDb = org.Hs.eg.db,
                                 ont = 'ALL',
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.05)
p21 <- dotplot(rna_table_enrich_res,label_format = 50)
saveplot(p21,'Fig5D',width = 4,height = 2)
atac_table <- data.table::fread('downstream_result/sample23/atac_snp/Beta_Type_2_diabetes.txt')
atac_table_enrich_res <- enrichGO(atac_table$symbol,
                                 keyType = 'SYMBOL',
                                 OrgDb = org.Hs.eg.db,
                                 ont = 'ALL',
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.05)
p22 <- dotplot(atac_table_enrich_res,label_format = 50)
saveplot(p22,'Fig5C',width = 4,height = 2)



map(total_rna_file,function(file_path){
    filename <- tools::file_path_sans_ext(file_path)
    atac_table <- data.table::fread(paste0(filename,'.txt'))
    atac_table_enrich_res <- enrichGO(atac_table$symbol,
                                      keyType = 'SYMBOL',
                                      OrgDb = org.Hs.eg.db,
                                      ont = 'ALL',
                                      pvalueCutoff = 0.05,
                                      qvalueCutoff = 0.05)
    p1 <- dotplot(atac_table_enrich_res)
    ggsave(filename = paste0(filename,'.png'),p1,width = 200,height = 100,units = 'mm',dpi = 300) 
},.progress = T)
########


#drct_table
total_samples <- list.files('downstream_result/')
for (Sample in total_samples) {
    drct_table <- data.table::fread(glue::glue('downstream_result/{Sample}/{Sample}_related_disease.xls'))
    sample_file <- map_dfr(1:nrow(drct_table),function(i){
        atac_searched_file <- list.files(glue::glue('downstream_result/{Sample}/atac_snp/'),pattern = gsub(' ','_',drct_table$Disease[i]),full.names = T) %>%
            stringr::str_subset('txt') %>%
            stringr::str_subset(gsub(' ','_',drct_table$`Cell type`[i])) 
        
        number_atac_files <- length(atac_searched_file)
        
        if (number_atac_files != 1) {
            atac_searched_file <- 'Unknown'
        }
        atac_df <- tibble(
            file = atac_searched_file,
            disease = paste0(drct_table$Disease[i]),
            ct = paste0(drct_table$`Cell type`[i]),
            type = 'ATAC',
            index = i
        )
        rna_searched_file <- list.files(glue::glue('downstream_result/{Sample}/rna_snp/'),pattern = gsub(' ','_',drct_table$Disease[i]),full.names = T) %>%
            stringr::str_subset('txt') %>%
            stringr::str_subset(gsub(' ','_',drct_table$`Cell type`[i])) 
        
        number_rna_files <- length(rna_searched_file)
        if (number_rna_files != 1) {
            rna_searched_file <- 'Unknown'
        }
        rna_df <- tibble(
            file = rna_searched_file,
            disease = paste0(drct_table$Disease[i]),
            ct = paste0(drct_table$`Cell type`[i]),
            type = 'RNA',
            index = i
        ) 
        df <- rbind(rna_df,atac_df)
        return(df)
    })
    drct_table <- drct_table[unique(sample_file[which(sample_file$file != 'Unknown'),]$index),]
    data.table::fwrite(drct_table,glue::glue('downstream_result/{Sample}/{Sample}_related_disease.xls'))
    sample_file <- map_dfr(1:nrow(drct_table),function(i){
        atac_searched_file <- list.files(glue::glue('downstream_result/{Sample}/atac_snp/'),pattern = gsub(' ','_',drct_table$Disease[i]),full.names = T) %>%
            stringr::str_subset('txt') %>%
            stringr::str_subset(gsub(' ','_',drct_table$`Cell type`[i])) 
        
        number_atac_files <- length(atac_searched_file)
        
        if (number_atac_files != 1) {
            atac_searched_file <- 'Unknown'
        }
        atac_df <- tibble(
            file = atac_searched_file,
            disease = paste0(drct_table$Disease[i]),
            ct = paste0(drct_table$`Cell type`[i]),
            type = 'ATAC',
            index = i
        )
        rna_searched_file <- list.files(glue::glue('downstream_result/{Sample}/rna_snp/'),pattern = gsub(' ','_',drct_table$Disease[i]),full.names = T) %>%
            stringr::str_subset('txt') %>%
            stringr::str_subset(gsub(' ','_',drct_table$`Cell type`[i])) 
        
        number_rna_files <- length(rna_searched_file)
        if (number_rna_files != 1) {
            rna_searched_file <- 'Unknown'
        }
        rna_df <- tibble(
            file = rna_searched_file,
            disease = paste0(drct_table$Disease[i]),
            ct = paste0(drct_table$`Cell type`[i]),
            type = 'RNA',
            index = i
        ) 
        df <- rbind(rna_df,atac_df)
        return(df)
    })
    if (length(which(sample_file$file == 'Unknown')) == 0) {
        data.table::fwrite(sample_file,glue::glue('downstream_result/{Sample}/{Sample}_snp_path.txt'),sep = '\t')
    }
    cat(Sample,'finished\n')
}

final_drct <- map_dfr(total_samples,function(Sample){
    drct_table <- data.table::fread(glue::glue('downstream_result/{Sample}/{Sample}_related_disease.xls')) %>%
        mutate(sample = Sample)
})


final_drct_snp_file <- map_dfr(total_samples,function(Sample){
    drct_table <- data.table::fread(glue::glue('downstream_result/{Sample}/{Sample}_snp_path.txt')) %>%
        mutate(sample = Sample)
})

final_drct_snp_rna_file <- final_drct_snp_file %>% filter(type == 'RNA')
final_drct_snp_atac_file <- final_drct_snp_file %>% filter(type == 'ATAC')


total_atac_file <- list.files('downstream_result',pattern = 'txt',recursive = T,full.names = T) %>%
    str_subset('atac_snp')
file.remove(total_atac_file[which(!(total_atac_file %in% final_drct_snp_atac_file$file))])

total_rna_file <- list.files('downstream_result',pattern = 'txt',recursive = T,full.names = T) %>%
    str_subset('rna_snp')
file.remove(total_rna_file[which(!(total_rna_file %in% final_drct_snp_rna_file$file))])



dar <- data.table::fread('downstream_result/sample23/sample23_scATAC-seq_95k_processed_DERs.txt.gz')
dar %>% filter(region == 'Promoter') %>% pull(gene)%>% unique() %>% length()
tf_act <- data.table::fread('downstream_result/sample23/sample23_tf_activity_top10.txt')

rna_snp <- data.table::fread('downstream_result/sample23/rna_snp/Beta_Type_2_diabetes.txt')
atac_snp <- data.table::fread('downstream_result/sample23/atac_snp/Beta_Type_2_diabetes.txt')

sample23_list <- list(DARs = dar,`TF activity` = tf_act,
                      `SNP overlapped peaks` = atac_snp, `SNP overlapped genes ` = rna_snp
                      )
writexl::write_xlsx(sample23_list,path = 'Supplementary Table 2.xlsx')

T2D <- data.table::fread('e://DRCTDb/data/LDSC_results/sample23/pvalues.tsv') %>%
    filter(V1 == 'Type_2_diabetes') %>%
    pivot_longer(-V1) %>%
    mutate(value = -log10(value)) %>%
    arrange(desc(value))
T2D$name[7] <- 'Immune'
T2D$name <- factor(T2D$name,levels = T2D$name )

p3 <- ggplot(T2D,aes(name,value)) +
    geom_col(aes(fill = name))+
    labs(x = NULL,y = '-log10(LDSC pvalue)',fill = NULL)  +
    geom_hline(yintercept = 3,linetype = 2) +
    scale_fill_manual(values = c("#982b2b","#db6968","#f47720", "#e5ce81","#0074b3","#88c4e8","#459943","#bdc3d2","#f0eedf") ) +
    theme_test()+
    guides(fill = 'none') +
    theme(axis.text.x = element_text(size = 8))
saveplot(p3,'Fig5A',width = 2.5,height = 1.5)


