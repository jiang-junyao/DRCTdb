library(Seurat)
library(tidyverse)
pbmc1 <- readRDS('../../data/scRNA-seq/Sample11/pbmc400k_cov19.rds')
pbmc2 <- readRDS('../../data/scRNA-seq/Sample11/pbmc750k_cov19.rds')

pbmc_healthy <- subset(pbmc1,COVID_status == 'Healthy')
saveRDS(pbmc_healthy,file = '../../data/scRNA-seq/Sample11/Sample11_processed_PBMC_healthy_173k.Rds')
