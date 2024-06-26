library(Seurat)
library(igraph)
library(GenomicRanges)
setwd('E:\\DRCTdb\\5.Downstream_analysis')
source('identify_region_motif.R')
source('overlap_gwas.R')
source('plot.R')
### path need to define
output_path = 'E:\\DRCTdb\\ignore\\downstream_result\\sample28\\'
rna_path = 'F:\\DRCTdb\\ignore\\scRNA-seq\\Sample28\\sample28_scRNA_46K_processed.Rds'
atac_path = 'E:\\DRCTdb\\ignore\\bed\\sample28/'
ldsc_path = "E:/DRCTdb/ignore/LDSC_results/sample28/pvalues.tsv"
#### db path

snp_path = 'E:\\public\\all_snp_info_gr.Rds'
disease_path = 'E:\\DRCTdb\\ignore\\LDSC_hg38\\summary_statistics\\Josh'
source('F:\\general_code\\run_cellchat.R')
source('F:\\general_code\\seurat_sample_ct.R')
### create output folder
dir.create(paste0(output_path))
dir.create(paste0(output_path,'grn_cor04'))
dir.create(paste0(output_path,'grn_cor02'))
dir.create(paste0(output_path,'rna_snp'))
dir.create(paste0(output_path,'atac_snp'))
dir.create(paste0(output_path,'ccc'))
###load data & define cell type use
rna = readRDS(rna_path)
ct = as.character(rna$cell_type)
ct = gsub(' ','_',ct)
ct = gsub('/','_',ct)
rna[['ct']] = ct
rna@active.ident = as.factor(rna$ct)

rna = seurat_sample_ct(rna,2000)
atac = dir(atac_path)
atac_ct = unlist(strsplit(atac,'.bed.gz'))
ct_use = atac_ct[c(1:5,7)]

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

ct_corres <- read.delim("E:/DRCTdb/ignore/LDSC_results/sample28/ct_corres.txt", header=FALSE)
### main calculation part
grn_list04 = list()
grn_list02 = list()
snp_list = list()
for (i in 1:length(ct_use)) {
  ### define significant disease
  rna_ct = ct_corres[ct_corres[,1]==ct_use[i],2]
  disease_use = ldsc_result[,1][ldsc_result[,ct_use[i]]<0.05]
  disease_use = unique(disease_use)
  disease_use = intersect(disease_use,names(disease_all))
  for (j in disease_use) {
      rna_use = subset(rna,ct==rna_ct)
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
### ldsc summary
ldsc_result = ldsc_result[!duplicated(ldsc_result[,1]),]
rownames(ldsc_result) = ldsc_result[,1]
ldsc_result = ldsc_result[,ct_use]
sig_ct_list = list()
sig_ct_vector = c()
for (i in 1:nrow(ldsc_result)) {
  sig_ct = colnames(ldsc_result)[ldsc_result[i,]<0.05]
  if (length(sig_ct)>0) {
    sig_ct_list[[rownames(ldsc_result)[i]]] = sig_ct
    sig_ct_merge = paste0(sig_ct,collapse = ';')
    sig_ct_vector = c(sig_ct_vector,sig_ct_merge)
  }
}
sig_ct_df = data.frame(names(sig_ct_list),sig_ct_vector)
colnames(sig_ct_df) = c('disease','related_cell_type')
write.table(sig_ct_df,paste0(output_path,'disease_related_celltypes.txt'),
            quote = F,sep = '\t',row.names = F)
### disease related ccc
source('F:\\general_code\\run_cellchat.R')
ccc_plot_list = list()
ccc_list = list()
for (i in 1:length(sig_ct_list)) {
  disease_name_use = names(sig_ct_list)[i]
  related_ct = sig_ct_list[[i]]
  rna_ct = ct_corres[ct_corres[,1] %in% related_ct,2]
  if (length(unique(rna_ct))>1) {
    rna_use = subset(rna,ct %in% rna_ct)
    ccc = run_cellchat(rna_use,rna_use@meta.data,group = 'ct',species = 'hs')
    groupSize <- as.numeric(table(ccc@idents))
    par(mfrow=c(1,1))
    p1 = netVisual_circle(ccc@net$weight, vertex.weight = groupSize,
                          weight.scale = T, label.edge= F,
                          title.name = "",vertex.label.cex = 0.7)
    print(disease_name_use)
    ### disease related ccc
    ccc_plot_list[[i]] = p1
    names(ccc_plot_list)[i] = names(sig_ct_list)[i]
    ccc_list[[i]] = ccc
    names(ccc_list)[i] = names(sig_ct_list)[i]
  }
}

### output ccc
for (i in 1:length(ccc_list)) {
  name_use = names(ccc_plot_list)[i]
  ccc = ccc_list[[name_use]]
  filenames = paste0(output_path,'ccc/',name_use,'.tiff')
  filenames2 = paste0(output_path,'ccc/',name_use,'.svg')
  if (!is.null(ccc_plot_list[[i]])) {
    svg(filename = filenames2, width = 4, height = 4)
    print(netVisual_circle(ccc@net$weight,
                           weight.scale = T, label.edge= F,
                           title.name = "",vertex.label.cex = 0.7))
    dev.off()
  }
}
saveRDS(ccc_list,paste0(output_path,'cellchat_ccc.rds'))

### output grn
for (i in 1:length(grn_list04)) {
  name_use = names(grn_list04)[i]
  out1 = grn_list04[[i]]
  write.table(out1,paste0(output_path,'grn_cor04/',name_use,'.txt')
              ,quote = F,sep = '\t',row.names = F)
  out1$type = ifelse(out1$value>0,'positive','negative')
  svg(filename = filenames2, width = 4, height = 4)
  print(ccc_plot_list[[i]])
  dev.off()
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
plot_heatmap_all()
plot_main(sample_use = 'sample27')

