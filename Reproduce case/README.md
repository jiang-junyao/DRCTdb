# A step-by-step guide for reproducing the data preparation process in DRCTDB

## Prepare seurat object of scRNA and scATAC

We using sample1 as example, first your should download all the rawdata form [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165838)

Run code in 01.Create_seurat_object.R

Will generate 2 seurat object (scRNA and scATAC) which contain peak-by-cell matrix and gene-by-cell matrix

Convert scATAC seurat object to h5ad for python analysis.

Many package can finish this. We using [dior](https://github.com/JiekaiLab/scDIOR) here.

After the conversion, we load it with scanpy to make sure its integrality(For all dataset, your can see the code in 1.Data_format_normalization/rds2h5ad.ipynb)




## Generate cell type specific peaks

Run code in 02.Generate_cell type_sepecific_peaks.R

> Note, your need to put the preprocess_functions.R (In this repository) in the same folder

This will generate a series bed file of cell type specific cCREs

## Run LDSC analysis
 
Just follow the example of [LDSC](https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses) or [HPC pipeline](https://github.com/kaizhang/LDSC)

your will obtained a table contain the p value matrix of each celltype against each disease (Check the format of pvalues.tsv)

>Note, if your have any problem in running LDSC, please contact us

## Obtain disease related cell types (DRCTs)

Run code in 03.Obtain_DRCTs.R

Will obtain a table contain DRCTs (Here we using p value <= 0.001 as threshold)

## Identity DERs and TF activity with scbasset

First 
>python run_scbasset.py sample1_scATAC-seq_80k_processed.h5ad

sample1_scATAC-seq_80k_processed.h5ad was geneated in step 01

This step will output a trained h5ad object,
a cell type differential open chromatin regions table, and a TF activity table

>Note, if your have any problem in running [scbasset](https://github.com/calico/scBasset), please contact us


Due to the scbasset TF activity is single cell level, we need to aggregate  TF activity in the same cell type (We select the top10 TFs in each cell type to make sure its correction)

Run 04.Get_TF_activity.R

your will obtained the 10 most active transcription factors in each cell type in Sample1_tf_activity.xls


