library(tidyverse)

Summary.Statistics <- read.csv("E:/Db/DRCTdb/test/Summary Statistics.txt", header=FALSE, row.names=NULL, sep=";")

df <- data.frame(
  trait = na.omit(str_squish(str_extract(Summary.Statistics$V1,'(?<==).*(?=,\\s)'))),
  file = na.omit(str_squish(str_extract(Summary.Statistics$V1,'(?<=trait_file =).*(?=\\s)')))
)
data.table::fwrite(df,'Summary.Statistics.xls',sep = '\t')

new_file <- data.frame(new_file = list.files('../../LDSC_hg38/summary_statistics/AlkesGroup/'))

data.table::fwrite(new_file,'Summary.Statistics2.xls',sep = '\t')

df <- readxl::read_excel('Summary.Statistics.xls',sheet = 'Sheet1')
df2 <- df %>% filter(!is.na(trait)) %>% select(1:2)

df2$file <- paste0('../summary_statistics/AlkesGroup/',df2$file)
data.table::fwrite(df2,'Summary.Statistics2.xls',sep = '\t')
df <- readxl::read_excel('Summary_Statistics.xlsx')
df2 <- df %>% distinct(.keep_all = T)

writexl::write_xlsx(df2,'Summary_Statistics2.xlsx')



df <- readxl::read_excel('Summary_Statistics.xlsx',sheet = 'Sheet1') %>% 
    select(2:1) %>% distinct(trait,.keep_all = T) %>% filter(!(trait %in% df2$trait))
df$new_file <- paste0('../summary_statistics/GWAS_Catalog/',df$new_file)



writexl::write_xlsx(df,'Summary_Statistics3.xlsx')

df2 <- readxl::read_excel('Summary_Statistics.xlsx',sheet = 'Sheet2')




