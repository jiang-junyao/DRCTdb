library(Seurat)
library(igraph)
library(GenomicRanges)
setwd('E:\\DRCTdb\\5.Downstream_analysis')
source('identify_region_motif.R')
source('overlap_gwas.R')
source('plot.R')
### path need to define
output_path = 'E:\\DRCTdb\\ignore\\downstream_result\\sample1\\'
rna_path = 'E:\\DRCTdb\\ignore\\scRNA-seq\\sample1\\sample1_heart_scRNA_35k_processed.Rds'
atac_path = 'E:\\DRCTdb\\ignore\\bed\\sample1/'

#### db path
ldsc_path = "E:/DRCTdb/ignore/LDSC_results/sample1/pvalues.tsv"
snp_path = 'E:\\public\\all_snp_info_gr.Rds'
disease_path = 'E:\\DRCTdb\\ignore\\LDSC_hg38\\summary_statistics\\Josh'

### create output folder
dir.create(paste0(output_path,'grn_cor04'))
dir.create(paste0(output_path,'grn_cor02'))
dir.create(paste0(output_path,'rna_snp'))
dir.create(paste0(output_path,'atac_snp'))

###load data & define cell type use
rna = readRDS(rna_path)
rna@active.ident=as.factor(rna$celltype)
atac = dir(atac_path)
atac_ct = unlist(strsplit(atac,'.bed.gz'))
ct_use = intersect(rna$celltype,atac_ct)
names(atac) = atac_ct
atac = atac[names(atac) %in% ct_use]
atac_list = list()
for (i in 1:length(atac)) {
  bed=read.delim(paste0(atac_path,atac[i]),header = F)
  peak = paste0(bed[,1],':',bed[,2],'-',bed[,3])
  atac_list[[names(atac)[i]]] = peak
}
ldsc_result <- read.delim(ldsc_path)
ct_use = intersect(colnames(ldsc_result),ct_use)
snp_all = readRDS(snp_path)
ldsc_result$X = gsub(' ','_',ldsc_result$X)

###get all disease
disease_all = dir(disease_path)
disease_all_name = as.data.frame(t(as.data.frame(strsplit(disease_all,'\\.'))))
names(disease_all) = disease_all_name[,4]


### main calculation part
grn_list04 = list()
grn_list02 = list()
snp_list = list()
for (i in 1:length(ct_use)) {
  ### define significant disease
  disease_use = ldsc_result[,1][ldsc_result[,ct_use[i]]<0.05]
  disease_use = unique(disease_use)
  disease_use = intersect(disease_use,names(disease_all))
  for (j in disease_use) {
      rna_use = subset(rna,celltype==ct_use[i])
      list1 = gwas_related_features(rna_use,atac_list[[i]],
                            disease_name=j,
                            snp_all = snp_all,zscore_thr = 1)
      
      grn04 = ct_grn_atac(list1[[2]][,1:3],
                        unique(list1[[1]]$symbol),
                        rna_use,cor_thr=0.4)
      grn02 = ct_grn_atac(list1[[2]][,1:3],
                        unique(list1[[1]]$symbol),
                        rna_use,cor_thr=0.2)
      grn_name = paste0(ct_use[i],'_',j)
      grn_list04[[grn_name]] = grn04
      grn_list02[[grn_name]] = grn02
      snp_list[[grn_name]] = list1
  }

}


### output grn
for (i in 1:length(grn_list04)) {
  name_use = names(grn_list04)[i]
  out1 = grn_list04[[i]]
  write.table(out1,paste0(output_path,'grn_cor04/',name_use,'.txt')
              ,quote = F,sep = '\t',row.names = F)
  out1$type = ifelse(out1$value>0,'positive','negative')
  plot_grn(out1)
}
for (i in 1:length(grn_list02)) {
  name_use = names(grn_list02)[i]
  out1 = grn_list02[[i]]
  write.table(out1,paste0(output_path,'grn_cor02/',name_use,'.txt')
              ,quote = F,sep = '\t',row.names = F)
}
### output snp
for (i in 1:length(snp_list)) {
  name_use = names(snp_list)[i]
  list_use = snp_list[[i]]
  rna_snp = list_use[[1]]
  atac_snp = list_use[[2]]
  write.table(rna_snp,paste0(output_path,'rna_snp/',name_use,'.txt'),
              quote = F,sep = '\t',row.names = F)
  write.table(atac_snp,paste0(output_path,'atac_snp/',name_use,'.txt'),
              quote = F,sep = '\t',row.names = F)
}
### visualize


########################
### infer grn from scRNA
########################
source('run_scenic.R')
run_scenic_main(rna_test)
rcistarget_out = readRDS('./int/2.1_tfModules_forMotifEnrichmet.Rds')
geneie3_out = readRDS('./int/1.4_GENIE3_linkList.Rds')
tf = t(as.data.frame(strsplit(names(rcistarget_out),'_')))[,1]
geneie3_out = geneie3_out[geneie3_out$TF %in% tf,]
geneie3_out = geneie3_out[geneie3_out$Target %in% unlist(rcistarget_out),]
