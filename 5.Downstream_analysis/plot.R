plot_grn_irena = function(data){
  library(igraph)
  library(IReNA)
  tf_use = unique(unlist(strsplit(Tranfac201803_Hs_MotifTFsF$TFs,';')))
  data_use = data[,c(1,2,5,3)]
  node_type = c(rep('TF',nrow(data_use)),
                ifelse(data_use$Var2 %in% tf_use,'TF','gene'))
  nodes = data.frame(c(data_use$Var1,data_use$Var2),node_type)
  colnames(nodes) <- c("name", "type")
  nodes <- nodes[!duplicated(nodes$name), ]
  colnames(data_use) <- c("from", "to", "type", "weight")
  g <- igraph::graph_from_data_frame(data_use, vertices = nodes, directed = TRUE)
  
  edge.color = c('#FDD1B0','#B3B3B3')
  ### define edge color
  edge.color2 <- c()
  for (i in data_use$type) {
    if (i == 'positive') {
      edge.color2 <- c(edge.color2,edge.color[1])
    } else if (i == 'negative') {
      edge.color2 <- c(edge.color2,edge.color[2])
    }
  }
  
  layout1 <- igraph::layout_on_grid(g)
  igraph::E(g)$arrow.size <- 0.8
  igraph::E(g)$arrow.width <- 0.5
  igraph::E(g)$label.color <- 'black'
  igraph::E(g)$color<- edge.color2
  igraph::E(g)$width<- 2.5
  igraph::V(g)$color <- ifelse(nodes$type=='TF','#67C7C1','#E56145')
  #igraph::V(g)$size<- vertex.size1
  igraph::V(g)$label.color <- 'black'
  igraph::V(g)$frame.color <- 'white'
  plot(g, layout = layout1, edge.curved = 0, vertex.label.cex =
         0.8, layout = layout1,
       vertex.shape='circle')
  if (length(unique(nodes$type))>1) {
    legend(x = 1.3, y = 1.3, levels(factor(igraph::V(g)$type)), pch = 21,
         col = "#777777", pt.bg = c('#E56145','#67C7C1'))
  }
  
}

plot_disease_heatmap <- function(pvalues){
  pvalues <- pvalues[!duplicated(pvalues[,1]),]
  rownames(pvalues) <- pvalues[,1]
  pvalues <- pvalues[,-1]
  p2 <- -log10(pvalues)

  Color1 <- c(rgb(102/255,46/255,115/255),rgb(31/255,153/255,139/255),rgb(251/255,232/255,48/255))

  pheatmap::pheatmap(as.matrix(p2), cluster_cols =T, cluster_rows =
                              T, color = colorRampPalette(Color1)(50),
                            border_color=rgb(200/255,200/255,200/255))
}
plot_heatmap_all <- function(ldsc_result_path='E:\\DRCTdb\\ignore\\LDSC_results/',
                             output_path = 'E:\\DRCTdb\\ignore\\downstream_result/'){
  dir1 = dir(output_path)
  for (i in dir1) {
    pvalue = read.delim(paste0(ldsc_result_path,i,'/pvalues.tsv'))
    out_path = paste0(output_path,i,'/ldsc_heatmap.tiff')
    out_path2 = paste0(output_path,i,'/ldsc_heatmap.svg')
    svg(filename = out_path2, width = 12, height = 10)
    print(plot_disease_heatmap(pvalue))
    dev.off()
  }
}
#plot_disease_heatmap(pvalues)

plot_enrich <- function(gene1){
  library(clusterProfiler)
  library(org.Hs.eg.db)
  entrez_id <- bitr(gene1, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)
  gene_universe <- names(org.Hs.eg.db)
  enrich_go <- enrichGO(gene = entrez_id$ENTREZID, 
                            keyType = 'ENTREZID', 
                            OrgDb = org.Hs.eg.db)
  p_go = dotplot(enrich_go)
  enrich_kegg <- enrichKEGG(gene = entrez_id$ENTREZID, 
                              organism = 'hsa',use_internal_data = T)
  p_kegg = dotplot(enrich_kegg)
  list1 = list(enrich_go,enrich_kegg)
  names(list1) = c('go','kegg')
  return(list1)
}
plot_main <- function(path_use,enrich=T,grn=T){
  dir1 = dir(path_use)
  for (i in dir1) {
    ### kegg go plot
    if (enrich) {
      
      library(ChIPseeker)
      library(TxDb.Hsapiens.UCSC.hg38.knownGene)
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
      library(org.Hs.eg.db)
      
      dir_rna = dir(paste0(path_use,i,'/rna_snp'))
      dir_atac = dir(paste0(path_use,i,'/atac_snp'))
      dir.create(paste0(path_use,i,'/rna_enrich_kegg'))
      dir.create(paste0(path_use,i,'/rna_enrich_go'))
      dir.create(paste0(path_use,i,'/atac_enrich_kegg'))
      dir.create(paste0(path_use,i,'/atac_enrich_go'))
      
      for (j in dir_rna) {
        snp_rna = read.delim(paste0(path_use,i,'/rna_snp/',j))
        enrich_plot = plot_enrich(snp_rna$symbol)
        name_disease = unlist(strsplit(j,'.txt'))
        kegg_path = paste0(path_use,i,'/rna_enrich_kegg/',name_disease,'_kegg.svg')
        go_path = paste0(path_use,i,'/rna_enrich_go/',name_disease,'_go.svg')
        
        if (dim(enrich_plot$kegg)[1]>1) {
          svg(filename = kegg_path, width = 8, height = 8)
          print(dotplot(enrich_plot$kegg))
          dev.off()
        }

        if (dim(enrich_plot$go)[1]>1) {
          svg(filename = go_path, width = 8, height = 8)
          print(dotplot(enrich_plot$go))
          dev.off()
        }
      }
      for (j in dir_atac) {
        snp_atac = read.delim(paste0(path_use,i,'/atac_snp/',j))
        snp_atac[,1] = paste0('chr',snp_atac[,1])
        peaks = GenomicRanges::GRanges(paste0(snp_atac[,1],':',
                                              snp_atac[,2],'-',
                                              snp_atac[,3]))
        peakAnno <- ChIPseeker::annotatePeak(peaks,
                                             tssRegion = c(-3000, 3000),
                                             TxDb = txdb, annoDb = 'org.Hs.eg.db')
        enrich_plot = plot_enrich(peakAnno@anno$SYMBOL)
        name_disease = unlist(strsplit(j,'.txt'))
        kegg_path = paste0(path_use,i,'/atac_enrich_kegg/',name_disease,'_kegg.tiff')
        go_path = paste0(path_use,i,'/atac_enrich_go/',name_disease,'_go.tiff')
        
        if (dim(enrich_plot$kegg)[1]>1) {
          svg(filename = kegg_path, width = 8, height = 8)
          print(dotplot(enrich_plot$kegg))
          dev.off()
        }
        
        if (dim(enrich_plot$go)[1]>1) {
          svg(filename = go_path, width = 8, height = 8)
          print(dotplot(enrich_plot$go))
          dev.off()
        }
      }
    }
    ### cell type specific grn
    if (grn) {
      print('plot grn!!')
      dir.create(paste0(path_use,i,'/grn_cor02_fig'))
      dir.create(paste0(path_use,i,'/grn_cor04_fig'))
      grn02_path = dir(paste0(path_use,i,'/grn_cor02'))
      grn04_path = dir(paste0(path_use,i,'/grn_cor04'))
      for (j in grn02_path) {
        grn_name = unlist(strsplit(j,'.txt'))
        data = read.table(paste0(path_use,i,'/grn_cor02/',j),header = T)
        data$type = ifelse(data$value>0,'positive','negative')
        if (nrow(data)>=1) {
          outpath = paste0(path_use,i,'/grn_cor02_fig/',grn_name,'.tiff')
          tiff(filename = outpath, width = 15000, height = 12000, units = "px", res = 1200, compression = "lzw")
          print(plot_grn_irena(data))
          dev.off()
        }
        
      }
      for (j in grn04_path) {
        grn_name = unlist(strsplit(j,'.txt'))
        data = read.table(paste0(path_use,i,'/grn_cor04/',j),header = T)
        data$type = ifelse(data$value>0,'positive','negative')
        if (nrow(data)>=1) {
          outpath = paste0(path_use,i,'/grn_cor04_fig/',grn_name,'.tiff')
          tiff(filename = outpath, width = 15000, height = 12000, units = "px", res = 1200, compression = "lzw")
          print(plot_grn_irena(data))
          dev.off()
        }
        
      }
    }

  }
}
