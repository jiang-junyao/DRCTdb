


output_path = '/storage/peiweikeLab/jiangjunyao/DRCTdb/ignore/test1/'
obj = readRDS('/storage/peiweikeLab/jiangjunyao/DRCTdb/ignore/test_obj.rds')
script_path = '/storage/peiweikeLab/jiangjunyao/DRCTdb/'
python3_path = '~/conda/bin/python'
original_config_path = '/storage/peiweikeLab/jiangjunyao/DRCTdb/ignore/LDSC_hg38/example/config.dhall'
ldsc_path = '/storage/peiweikeLab/jiangjunyao/DRCTdb/ignore/LDSC_hg38/ldsc-Ubuntu-x86_64'
ncore = 10

###1.marker peak
source(paste0(script_path,'2.Data preprocessing/preprocess.R'))
DBRs_gr <- find_DBRs(obj)
peak_bed = as.data.frame(DBRs_gr)
cluster_use =levels(obj@active.ident)
cluster_use_path = paste0(output_path,'cluster_use.txt')
write.table(cluster_use,paste0(),quote = F,
            sep = '\t',col.names = F,row.names = F)
for (i in cluster_use) {
  peak_bed_cluster = peak_bed[peak_bed$cluster %in% i,]
  write.table(peak_bed_cluster,paste0(output_path,i,'.bed'))
}
###2.ldsc
config_out_path = paste0(output_path,'config_use.dhall')
system(paste0(python3_path,' ',script_path,
              '3.Identify_disease_related_cell_type/change_config.py ',
              cluster_use_path,' ',original_config_path,' ',config_out_path))
system(paste0(ldsc_path,' run --config ',config_out_path,' -n ',ncore))

###3.scRNA marker

###4.downstream