library(readxl)

# 数据加载
coretable <- readxl::read_excel('scdb_core.xlsx', sheet = 'Sheet3')
download_table <- readxl::read_excel('downlaod_link.xlsx')

# 数据预处理
download_table$scRNA <- download_table$`Processed scRNA`
download_table$`Processed scRNA` <- 'Download data'
download_table$scATAC <- download_table$`Processed scATAC`
download_table$`Processed scATAC` <- 'Download data'
download_table$ldsc <- download_table$`LDSC results`
download_table$`LDSC results` <- 'Download data'

# 创建下载链接
download_table$`Processed scRNA` <- sprintf('<a href="%s" target="_blank">%s</a>', download_table$scRNA, download_table$`Processed scRNA`)
download_table$`Processed scATAC` <- sprintf('<a href="%s" target="_blank">%s</a>', download_table$scATAC, download_table$`Processed scATAC`)
download_table$`LDSC results` <- sprintf('<a href="%s" target="_blank">%s</a>', download_table$ldsc, download_table$`LDSC results`)
