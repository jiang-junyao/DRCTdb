library(tidyverse)
df <- data.table::fread('pvalues.tsv') %>% as_tibble() %>% 
    pivot_longer(cols = -V1,names_to = 'cell_type',values_to = 'pvalue') %>% 
    mutate(neglog10p = -log10(pvalue)) %>% 
    arrange(desc(neglog10p)) %>%
    filter(V1 != 'Atrial_fibrillation',abs(neglog10p) >= 3)  %>% 
    dplyr::select(1:3) 
colnames(df) <- c('Disease', 'Cell type','p value')
data.table::fwrite(df, file = 'Sample1_related_disease.xls',sep = '\t')
