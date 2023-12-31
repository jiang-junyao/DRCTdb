###load data
library(Seurat)
test = Read10X('E:\\public\\public_data\\10X\\multiomics\\PBMC_from_a Healthy_Donor_No_Cell_Sorting_(3k)\\pbmc_unsorted_3k_filtered_feature_bc_matrix\\filtered_feature_bc_matrix')
peaks = test$Peaks
peaks = rownames(peaks)[1:10000]
peaks = strsplit(peaks,':')
peaks = as.data.frame(t(as.data.frame(peaks)))
region = strsplit(peaks$V2,'-')
region = as.data.frame(t(as.data.frame(region)))
peaks = cbind(peaks[,1],region)
rna_test = test$`Gene Expression`
rna_test = rna_test[1:1000,1:1000]
rna_test = as.matrix(rna_test)
colnames(peaks) = c('V1','V2','V3')

########################
###identify tf-target from scatac
########################
library(IReNA)
source('identify_region_motif.R')
### peak annotation
peak_gr = peak_anno(peaks)
### find peak related motif
PWM = readRDS('F:\\public\\Transfac_PWMatrixList.rds')
motif1 = Tranfac201803_Hs_MotifTFsF
gene.use = rownames(rna_test)
peak_motif_gr = identify_region_tfs(peak_gr,gene.use,PWM,motif1)
### overlap motif and peak
atac_out = overlap_peak_motif(peak_gr,peak_motif_gr,motif1)
tf_target = make_tf_target(atac_out)
########################
### infer grn from scRNA
########################
source('run_scenic.R')
run_scenic_main(rna_test)
rcistarget_out = readRDS('./int/2.1_tfModules_forMotifEnrichmet.Rds')
geneie3_out = readRDS('./int/1.4_GENIE3_linkList.Rds')
tf = t(as.data.frame(strsplit(names(rcistarget_out),'_')))[,1]
geneie3_out = geneie3_out[geneie3_out$TF %in% tf,]
geneie3_out = geneie3_out[geneie3_out$Target %in% unlist(rcistarget_out),]

########################
### integration
########################
geneie3_out$idx = paste0(geneie3_out$TF,'-',geneie3_out$Target)
geneie3_out = geneie3_out[geneie3_out$idx %in% tf_target,]





library(tidyverse)

transform_value <- function(vector){
    num_values <- length(vector)
    top_10_percent_index <- ceiling(0.1 * num_values)
    top_50_percent_index <- ceiling(0.5 * num_values)
    
    sorted_indices <- order(vector, decreasing = TRUE)
    
    vector[sorted_indices[1:top_10_percent_index]] <- 20
    
    
    vector[sorted_indices[(top_10_percent_index + 1):top_50_percent_index]] <- 10
    
    vector[sorted_indices[(top_50_percent_index + 1):num_values]] <- 5
    return(vector)
}

get_GRN_NODE <- function(df){
    node_info <- tibble(
        gene = unique(c(df$Var1,df$Var2))
    )
    node_info$group <- if_else(node_info$gene %in% df$Var1,'TF','Not TF')
    node_info$group2 <- c(c(1:sum(node_info$group == 'TF')),
                          c(rep(sum(node_info$group == 'TF') + 1,sum(node_info$group != 'TF'))))
    node_info$size <- map_int(node_info$gene,function(x){
        number <- sum(str_count(df$Var1,x))
    })
    node_info$size <- transform_value(node_info$size)
    df <- as.data.frame(df)
    node_info <- as.data.frame(node_info)
    gene_vec <- 0:(length(node_info$gene)-1)
    names(gene_vec) <- node_info$gene
    df$source <- map_int(df$Var1,function(x){
        gene_vec[x]
    })
    df$target <- map_int(df$Var2,function(x){
        gene_vec[x]
    })
    
    grn_list <- list(
        link = df,
        nodes = node_info
    )
    return(grn_list)
}

df <- data.table::fread('./downstream_result/sample1/grn_cor02/aCM_Atrial_fibrillation.txt')
df <- data.table::fread('./downstream_result/sample1/grn_cor02/aCM_Myocardial_fractal_dimension_(slices_1-9).txt')
aCM_Atrial_fibrillation <- get_GRN_NODE(df)


grn_file <- list.files('.',pattern = 'txt',recursive = T,full.names = T) %>%
    str_subset('grn_cor02')


map(grn_file,function(i){
    filename <- tools::file_path_sans_ext(i)
    rds_path <- paste0(filename,'.Rds')
    grn_table <- data.table::fread(i)
    grn_table <- grn_table %>% arrange(desc(value))
    if (nrow(grn_table) > 1000) {
        grn_table <- grn_table[1:1000,]
    }
    if (nrow(grn_table) > 1) {
        grn_lsit <- get_GRN_NODE(grn_table)
        saveRDS(grn_lsit,rds_path)
    }
},.progress = T)


test <- readRDS('downstream_result/sample10/grn_cor02/aDC_Area_under_the_curve_of_insulin_levels.Rds')

# library(networkD3)
# forceNetwork(Links = test$link, Nodes = test$nodes,
#              Source = "source", Target = "target",fontSize = 10,
#              Value = "value", NodeID = "gene",zoom = TRUE,arrows =TRUE,
#              Group = "group2", opacity = 0.8)