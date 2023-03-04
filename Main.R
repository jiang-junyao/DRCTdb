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
