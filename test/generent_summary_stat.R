library(tidyverse)

Summary.Statistics <- read.csv("E:/Db/DRCTdb/test/Summary Statistics.txt", header=FALSE, row.names=NULL, sep=";")

df <- data.frame(
  trait = na.omit(str_squish(str_extract(Summary.Statistics$V1,'(?<==).*(?=,\\s)'))),
  file = na.omit(str_squish(str_extract(Summary.Statistics$V1,'(?<=trait_file =).*(?=\\s)')))
)
data.table::fwrite(df,'Summary.Statistics.xls',sep = '\t')

new_file <- data.frame(new_file = list.files('../../LDSC_hg38/summary_statistics/AlkesGroup/'))

data.table::fwrite(new_file,'Summary.Statistics2.xls',sep = '\t')
