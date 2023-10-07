library(Seurat)
library(GenomicRanges)
### part1
merged_matrix <- ReadMtx('../ignore/GSE165837_CARE_ATAC_merged_matrix.mtx.gz',
                         cells = '../ignore/GSE165837_CARE_ATAC_merged_barcodes.txt.gz',
                         features = '../ignore/GSE165837_CARE_ATAC_merged_features.txt.gz',
                         feature.column = 1)
label <- read.delim("E:/DRCTdb/ignore/GSE165837_CARE_ATAC_merged_umap_cluster_labels.tsv.gz", row.names=1)
obj <- CreateSeuratObject(merged_matrix,meta.data = label)
save(merged_matrix,file = '../data/scATAC-seq/sample1/merged_matrix.Rds')
### part2
source('2.Data preprocessing/preprocess.R')
DBRs_gr <- find_DBRs(obj)
saveRDS(DBRs_gr,'../ignore/cardiac_DBRs_wilcox_gr.rds')

### part3
### Note: please formatting these eqtl data to get official inputs 
### (e.g. chromosome name)
### This part is just the example
### Due to the data, i have not test the performance
source('3.Identify_disease_related_cell_type/cal_ct_score.R')
eqtl <- read.delim("../ignore/GTEx_v8_finemapping_CAVIAR/CAVIAR_Results_v8_GTEx_LD_HighConfidentVariants.gz")
eqtl <- eqtl[eqtl$TISSUE %in% 'Heart_Atrial_Appendage',]
eqtl$CHROM <- paste0('chr',eqtl$CHROM)
eqtl_gr <- GRanges(paste0(eqtl$CHROM,':',eqtl$POS,'-',eqtl$POS))
overlapregion <- findOverlaps(DBRs_gr,eqtl_gr)
overlap_DBRs <- DBRs_gr[overlapregion@from]
ct_score <- cal_ct_score(overlap_DBRs,DBRs_gr)



##Basic statistics
library(tidyverse)
sample_tissue <- readxl::read_excel('../data/sample_tissue.xlsx',col_names = F)
colnames(sample_tissue) <- c('dataset','tissue')
sample_tissue <- separate_rows(sample_tissue,tissue,sep = ',')
writexl::write_xlsx(sample_tissue,'../data/sample_tissue2.xlsx')

sample1 <- readRDS('../data/Rds/sample1_scATAC-seq_80k_processed.Rds')
sample3 <- readRDS('../data/Rds/sample3_scATAC-seq_756k_processed.Rds')

sample3_df <- sample3$cell_type |> table()  |> as.data.frame()
sample3_meta <- data.table::fread('../data/scATAC-seq/Sample3/Cell_metadata.tsv.gz')

sample3_tissue <- sample_tissue %>% filter(dataset == 'Sample3') %>% pull(tissue)
map_int(sample3_tissue,function(x){
    sample3_df %>% filter(str_detect(sample3_df$Var1,tolower(x))) %>% pull(Freq) %>% sum()
}) %>% setNames(sample3_tissue)

sample3_df2 <- sample3_meta$tissue |> table() |> as.data.frame()
sample3_df2 %>% filter(str_detect(Var1,tolower('Breast'))) %>% pull(Freq) %>% sum()
sample4 <- readRDS('../data/Rds/sample4_scATAC-seq_30k_processed.Rds')

sample16 <- readRDS('../data/Rds/sample16_Bone_marrow_ATAC_Healthy_35k_processed.Rds')


##
sample5_meta <- data.table::fread('../data/scATAC-seq/Sample5/filtered.cell_metadata.for_website.txt.gz')
sample26 <- readRDS('../data/Rds/Sample26_scATAC-seq_229k_processed.Rds')
sample26_df <- sample26$Sample |> table()  |> as.data.frame()
sample26_df %>% filter(str_detect(Var1,'BY')) %>% pull(Freq) %>% sum()



sample_tissue <- readxl::read_excel('../data/scdb_core.xlsx',sheet = 'Sheet2')

test <- sample_tissue %>% group_by(dataset) %>% summarise(tissue = paste0(tissue,collapse = ';')) %>% 
    mutate(num = str_extract(dataset,'\\d+'))
    arrange(dataset, .by_group = TRUE)
