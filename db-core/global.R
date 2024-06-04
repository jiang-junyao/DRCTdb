library(readxl)
library(enrichR)



# 数据加载
coretable <- readxl::read_excel('scdb_core.xlsx', sheet = 'Sheet3')
download_table <- readxl::read_excel('downlaod_link.xlsx')

# 数据预处理
download_table$scRNA <- download_table$`Processed scRNA`
download_table$`Processed scRNA` <- 'Download data'
download_table$scATAC <- download_table$`Processed scATAC`
download_table$`Processed scATAC` <- 'Download data'
#download_table$ldsc <- download_table$`LDSC results`
#download_table$`LDSC results` <- 'Download data'

enrich_GO <- function(df){
    dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023")
    enriched <- enrichR::enrichr(df[['SYMBOL']], dbs)
    return(enriched)
}

# 创建下载链接
download_table$`Processed scRNA` <- sprintf('<a href="%s" target="_blank">%s</a>', download_table$scRNA, download_table$`Processed scRNA`)
download_table$`Processed scATAC` <- sprintf('<a href="%s" target="_blank">%s</a>', download_table$scATAC, download_table$`Processed scATAC`)

#download_table$`LDSC results` <- sprintf('<a href="%s" target="_blank">%s</a>', download_table$ldsc, download_table$`LDSC results`)


#library(org.Hs.eg.db)
#AnnotationDbi::saveDb(org.Hs.eg.db, "org_HS.sqlite")
#org.Hs.eg.db <- AnnotationDbi::loadDb('org_HS.sqlite')
#PDPN





