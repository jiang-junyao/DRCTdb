###load data
library(Seurat)
test = Read10X('E:\\public\\public_data\\10X\\multiomics\\PBMC_from_a Healthy_Donor_No_Cell_Sorting_(3k)\\pbmc_unsorted_3k_filtered_feature_bc_matrix\\filtered_feature_bc_matrix')
peaks = test$Peaks
peaks = rownames(peaks)[1:10000]
peaks = strsplit(peaks,':')
peaks = as.data.frame(t(as.data.frame(peaks)))
region = strsplit(peaks$V2,'-')
region = as.data.frame(t(as.data.frame(region)))
peaks = cbind(peaks[,1],region)
rna_test = test$`Gene Expression`
rna_test = rna_test[1:1000,1:1000]
rna_test = as.matrix(rna_test)
colnames(peaks) = c('V1','V2','V3')

########################
###identify tf-target from scatac
########################
library(IReNA)
source('identify_region_motif.R')
### peak annotation
peak_gr = peak_anno(peaks)
### find peak related motif
PWM = readRDS('F:\\public\\Transfac_PWMatrixList.rds')
motif1 = Tranfac201803_Hs_MotifTFsF
gene.use = rownames(rna_test)
peak_motif_gr = identify_region_tfs(peak_gr,gene.use,PWM,motif1)
### overlap motif and peak
atac_out = overlap_peak_motif(peak_gr,peak_motif_gr,motif1)
tf_target = make_tf_target(atac_out)
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

########################
### integration
########################
geneie3_out$idx = paste0(geneie3_out$TF,'-',geneie3_out$Target)
geneie3_out = geneie3_out[geneie3_out$idx %in% tf_target,]
