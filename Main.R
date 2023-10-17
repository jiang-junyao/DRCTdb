library(Seurat)
library(GenomicRanges)
### part1
merged_matrix <- ReadMtx('../ignore/GSE165837_CARE_ATAC_merged_matrix.mtx.gz',
                         cells = '../ignore/GSE165837_CARE_ATAC_merged_barcodes.txt.gz',
                         features = '../ignore/GSE165837_CARE_ATAC_merged_features.txt.gz',
                         feature.column = 1)
label <- read.delim("E:/DRCTdb/ignore/GSE165837_CARE_ATAC_merged_umap_cluster_labels.tsv.gz", row.names=1)
obj <- CreateSeuratObject(merged_matrix,meta.data = label)
save(merged_matrix,file = '../data/scATAC-seq/sample1/merged_matrix.Rds')
### part2
source('2.Data preprocessing/preprocess.R')
DBRs_gr <- find_DBRs(obj)
saveRDS(DBRs_gr,'../ignore/cardiac_DBRs_wilcox_gr.rds')

### part3
### Note: please formatting these eqtl data to get official inputs 
### (e.g. chromosome name)
### This part is just the example
### Due to the data, i have not test the performance
source('3.Identify_disease_related_cell_type/cal_ct_score.R')
eqtl <- read.delim("../ignore/GTEx_v8_finemapping_CAVIAR/CAVIAR_Results_v8_GTEx_LD_HighConfidentVariants.gz")
eqtl <- eqtl[eqtl$TISSUE %in% 'Heart_Atrial_Appendage',]
eqtl$CHROM <- paste0('chr',eqtl$CHROM)
eqtl_gr <- GRanges(paste0(eqtl$CHROM,':',eqtl$POS,'-',eqtl$POS))
overlapregion <- findOverlaps(DBRs_gr,eqtl_gr)
overlap_DBRs <- DBRs_gr[overlapregion@from]
ct_score <- cal_ct_score(overlap_DBRs,DBRs_gr)






##------
scdb_core <- readxl::read_xlsx('../../data/scdb_core.xlsx')
scdb_core_total <- scdb_core %>% separate_rows(total_disease,sep = ';') %>% separate_rows(Tissue,sep = ';')
writexl::write_xlsx(scdb_core_total,'../../data/scdb_core_total.xlsx')
##Basic statistics
library(tidyverse)
sample_tissue <- readxl::read_excel('../data/sample_tissue.xlsx',col_names = F)
colnames(sample_tissue) <- c('dataset','tissue')
sample_tissue <- separate_rows(sample_tissue,tissue,sep = ',')
writexl::write_xlsx(sample_tissue,'../data/sample_tissue2.xlsx')

sample1 <- readRDS('../data/Rds/sample1_scATAC-seq_80k_processed.Rds')
sample3 <- readRDS('../data/Rds/sample3_scATAC-seq_756k_processed.Rds')

sample3_df <- sample3$cell_type |> table()  |> as.data.frame()
sample3_meta <- data.table::fread('../data/scATAC-seq/Sample3/Cell_metadata.tsv.gz')

sample3_tissue <- sample_tissue %>% filter(dataset == 'Sample3') %>% pull(tissue)
map_int(sample3_tissue,function(x){
    sample3_df %>% filter(str_detect(sample3_df$Var1,tolower(x))) %>% pull(Freq) %>% sum()
}) %>% setNames(sample3_tissue)

sample3_df2 <- sample3_meta$tissue |> table() |> as.data.frame()
sample3_df2 %>% filter(str_detect(Var1,tolower('Breast'))) %>% pull(Freq) %>% sum()
sample4 <- readRDS('../data/Rds/sample4_scATAC-seq_30k_processed.Rds')

sample16 <- readRDS('../data/Rds/sample16_Bone_marrow_ATAC_Healthy_35k_processed.Rds')


##
sample5_meta <- data.table::fread('../data/scATAC-seq/Sample5/filtered.cell_metadata.for_website.txt.gz')
sample26 <- readRDS('../data/Rds/Sample26_scATAC-seq_229k_processed.Rds')
sample26_df <- sample26$Sample |> table()  |> as.data.frame()
sample26_df %>% filter(str_detect(Var1,'BY')) %>% pull(Freq) %>% sum()



sample_tissue <- readxl::read_excel('../data/scdb_core.xlsx',sheet = 'Sheet2')

test <- sample_tissue %>% group_by(dataset) %>% summarise(tissue = paste0(tissue,collapse = ';')) %>% 
    mutate(num = as.numeric(str_extract(dataset,'\\d+'))) %>% 
    arrange(num)

sample_tissue2 <- sample_tissue %>% group_by(tissue) %>%  summarise(cell_num = sum(Cell_num))


doughnut <- function (x, labels = names(x), edges = 200, outer.radius = 0.8,
                      inner.radius=0.6, clockwise = FALSE,
                      init.angle = if (clockwise) 90 else 0, density = NULL,
                      angle = 45, col = NULL, border = FALSE, lty = NULL,
                      main = NULL, ...)
{
    if (!is.numeric(x) || any(is.na(x) | x < 0))
        stop("'x' values must be positive.")
    if (is.null(labels))
        labels <- as.character(seq_along(x))
    else labels <- as.graphicsAnnot(labels)
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    plot.new()
    pin <- par("pin")
    xlim <- ylim <- c(-1, 1)
    if (pin[1L] > pin[2L])
        xlim <- (pin[1L]/pin[2L]) * xlim
    else ylim <- (pin[2L]/pin[1L]) * ylim
    plot.window(xlim, ylim, "", asp = 1)
    if (is.null(col))
        col <- if (is.null(density))
            palette()
    else par("fg")
    col <- rep(col, length.out = nx)
    border <- rep(border, length.out = nx)
    lty <- rep(lty, length.out = nx)
    angle <- rep(angle, length.out = nx)
    density <- rep(density, length.out = nx)
    twopi <- if (clockwise)
        -2 * pi
    else 2 * pi
    t2xy <- function(t, radius) {
        t2p <- twopi * t + init.angle * pi/180
        list(x = radius * cos(t2p),
             y = radius * sin(t2p))
    }
    for (i in 1L:nx) {
        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(x[i], x[i + 1], length.out = n),
                  outer.radius)
        polygon(c(P$x, 0), c(P$y, 0), density = density[i],
                angle = angle[i], border = border[i],
                col = col[i], lty = lty[i])
        Pout <- t2xy(mean(x[i + 0:1]), outer.radius)
        lab <- as.character(labels[i])
        if (!is.na(lab) && nzchar(lab)) {
            lines(c(1, 1.05) * Pout$x, c(1, 1.05) * Pout$y)
            text(1.1 * Pout$x, 1.1 * Pout$y, labels[i],
                 xpd = TRUE, adj = ifelse(Pout$x < 0, 1, 0),
                 ...)
        }      
        Pin <- t2xy(seq.int(0, 1, length.out = n*nx),
                    inner.radius)
        polygon(Pin$x, Pin$y, density = density[i],
                angle = angle[i], border = border[i],
                col = "white", lty = lty[i])
    }
    
    title(main = main, ...)
    invisible(NULL)
}

clustcol<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
cell_num <- shuffle(sample_tissue2$cell_num)
names(cell_num) <- sample_tissue2$tissue


pdf('../data/tissue_cell_num.pdf',width = 8,height = 8)
doughnut(sample(cell_num),col = clustcol)
dev.off()

for (i in list.files('../../data/downstream_result/')) {
    raw_path <- paste0('../../data/downstream_result/',i,'/cellchat_ccc.rds')
    new_path <- paste0('../../data/CCC_obj/',i,'_ccc.rds')
    file.copy(raw_path,new_path,overwrite = T)
    file.remove(raw_path)
    cat(i,'finished\n')
}
for (i in list.files('../../data/DERs/')){
    raw_path <- paste0('../../data/DERs/',i)
    new_path <- paste0('../../data/downstream_result/',str_extract(i,'sample\\d+'),'/',i)
    file.copy(raw_path,new_path,overwrite = T)
}

for (i in list.files('../../data/LDSC_results/',pattern = 'zip')){
    raw_path <- paste0('../../data/LDSC_results/',i)
    new_path <- paste0('../../data/downstream_result/',str_extract(i,'sample\\d+'),'/',i)
    file.copy(raw_path,new_path,overwrite = T)
    file.remove(raw_path)
    cat(i,'finished\n')
}

##text to json

library(jsonlite)


table2json <- function(x){
    filename <- tools::file_path_sans_ext(basename(x))
    output_name <- paste0(dirname(x),'/',filename,'.json')
    if (str_detect(x,'xlsx')) {
        df <- readxl::read_excel(x)
    }else{
        df <- data.table::fread(x)
    }
    json_data <- jsonlite::toJSON(df)
    jsonlite::write_json(json_data,output_name)
    return(json_data)
}


#test <-  table2json('../../data/scdb_core.xlsx')


all_downstream_txt <- list.files('../data/downstream_result/',pattern = 'txt|xls',recursive = T,full.names = T)


res <- map(all_downstream_txt,table2json,.progress = T)
