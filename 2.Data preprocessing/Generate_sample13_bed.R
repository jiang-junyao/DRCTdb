library(Signac)
library(Seurat)
library(tidyverse)
library(Matrix)
library(glue)
source('preprocess_functions.R')
sample13_ATAC <- readRDS('../../data/scATAC-seq/Sample13/sample13_Bone_marrow_ATAC_4k_processed.Rds.rds')
metadata <- sample13_ATAC@meta.data
