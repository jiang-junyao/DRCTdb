
gwas_related_features <- function(rna_use,atac_use,snp_all,
                                  gwas_path='E:\\DRCTdb\\ignore\\LDSC_hg38/',
                                  disease_name,
                                  zscore_thr = 2){
  rna_features = rownames(extract_expressed_features(rna_use,cells_quantile = 0.1))
  atac_features = rownames(extract_expressed_features(atac_use,cells_quantile =0.025))
  atac_features = gsub('chr','',atac_features)
  peak_gr = GenomicRanges::GRanges(atac_features)
  gene_tss = readRDS('E:\\public/hg38_gene_tss.rds')
  gene_tss_use = gene_tss[gene_tss$symbol %in% rna_features,]
  gene_tss_region = get_tss(gene_tss_use)
  gene_tss_region[,1] = gsub('chr','',gene_tss_region[,1])
  gene_gr = GenomicRanges::GRanges(paste0(gene_tss_region[,1],':',
                                          gene_tss_region[,2],'-',
                                          gene_tss_region[,3]))
  gene_gr$symbol = gene_tss_region[,4]
  gene_overlap = overlap_gwas(gene_gr,gwas_path = gwas_path
                              ,snp_all = snp_all,disease_name = disease_name,
                              zscore_thr = zscore_thr)
  
  peak_overlap = overlap_gwas(peak_gr,gwas_path = gwas_path
                              ,snp_all = snp_all,disease_name = disease_name,
                              zscore_thr = zscore_thr)
  return(list(gene_overlap,peak_overlap))
}

extract_expressed_features <- function(obj,cells_quantile = 0.05){
  matrix_use = GetAssayData(obj)
  if (cells_quantile==0) {
    TfExp <- matrix_use[rowSums(as.matrix(matrix_use))>0,]
  }else{
    quantile_exp <- ncol(as.matrix(matrix_use))/(1/cells_quantile)
    TfExp <- matrix_use[ncol(as.matrix(matrix_use))-rowSums(as.matrix(matrix_use==0))>quantile_exp,]}
  return(TfExp)
}

overlap_gwas <- function(peak_gr,gwas_path = 'E:\\DRCTdb\\ignore\\LDSC_hg38/',
                         snp_all,
                        disease_name,zscore_thr = 2){
  ###load data
  josh_path = paste0(gwas_path,'summary_statistics/Josh/')
  disease_all = dir(josh_path)
  disease_all_name = as.data.frame(t(as.data.frame(strsplit(disease_all,'\\.'))))
  names(disease_all) = disease_all_name[,4]
  
  gwas_file = disease_all[disease_name]
  print(gwas_file)
  gwas_file = read.delim(paste0(josh_path,gwas_file))
  
  gwas_file = gwas_file[!is.na(gwas_file[,5]),]
  gwas_file = gwas_file[gwas_file$Z > zscore_thr | gwas_file$Z < (-zscore_thr),]
  snp_use = intersect(gwas_file$SNP,snp_all$snp)
  gwas_file = gwas_file[gwas_file$SNP %in% snp_use,]
  snp_use = snp_all[match(gwas_file$SNP,snp_all$snp)]
  overlapped_peak = findOverlaps(peak_gr,snp_use)
  peak_out = cbind(as.data.frame(peak_gr[overlapped_peak@from,]),
                   as.data.frame(snp_use[overlapped_peak@to,]))
  return(peak_out)
}

get_tss <- function(final,upstream_length=1000,downstream_length=500){
  start1 <- c()
  end1 <- c()
  for (i in 1:nrow(final)) {
    if (final[i,4]=='-') {
      start1 <- c(start1,as.numeric(final[i,3])-upstream_length)
      end1 <- c(end1,as.numeric(final[i,3])+downstream_length)
    }else{
      start1 <- c(start1,as.numeric(final[i,2])-upstream_length)
      end1 <- c(end1,as.numeric(final[i,2])+downstream_length)
    }
  }
  output1 = data.frame(final[,1],start1,end1,final[,6])
  colnames(output1) = c('chr','start','end','gene')
  return(output1)
}
