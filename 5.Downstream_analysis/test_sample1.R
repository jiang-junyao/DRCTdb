###load data
library(Seurat)
library(igraph)
setwd('E:\\DRCTdb\\5.Downstream_analysis')
source('identify_region_motif.R')
source('overlap_gwas.R')
rna = readRDS('E:\\DRCTdb\\ignore\\sample1_heart_scRNA_35k_processed.Rds')
atac = readRDS('E:\\DRCTdb\\ignore\\sample1_seurat_obj.Rds')
#ct_use = intersect(rna$celltype,atac$cell_type)
ct_use = c('aCM','vCM','EC')
disease_use = c('Atrial_fibrillation','Atrial_fibrillation','Varicose_veins')
snp_all = readRDS('E:\\public\\all_snp_info_gr.Rds')
peak_name = t(as.data.frame(strsplit(rownames(atac),'-')))
peak_name = paste0(peak_name[,1],':',peak_name[,2],'-',peak_name[,3])
grn_list = list()
for (i in 1:length(ct_use)) {
  
  rna_use = subset(rna,celltype==ct_use[i])
  atac_use = subset(atac,cell_type==ct_use[i])
  rownames(atac_use@assays$peaks@counts) = peak_name
  rownames(atac_use@assays$peaks@data) = peak_name
  list1 = gwas_related_features(rna_use,atac_use,
                        disease_name=disease_use[3],
                        snp_all = snp_all,zscore_thr = 1)
  
  grn = ct_grn_atac(list1[[2]][,1:3],
                    unique(list1[[1]]$symbol),
                    rna_use,cor_thr=0.2)
  grn_list[[i]] = grn
}
### visualize
for (i in grn_list) {
  g = graph_from_data_frame(i[,1:3])
  igraph::V(g)$size = 20
  igraph::E(g)$arrow.size = 0.5
  igraph::E(g)$arrow.width = 0.5
  layout1 <- layout_on_grid(g)
  plot(g,edge.curved = 0,layout = layout1,vertex.label.cex =0.8)
}

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
