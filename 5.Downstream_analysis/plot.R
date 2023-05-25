setwd('E:\\DRCTdb\\5.Downstream_analysis')
data = read.table('../ignore/downstream_result/sample1/grn_cor04/aCM_Atrial_fibrillation.txt',header = T)
data$type = ifelse(data$value>0,'positive','negative')
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
  legend(x = 1.3, y = 1.3, levels(factor(igraph::V(g)$type)), pch = 21,
         col = "#777777", pt.bg = c('#E56145','#67C7C1'))
}
#plot_grn_irena(data)



plot_grn <- function(data){
  library(ggraph)
  ggraph(data,layout = 'fr') + 
    geom_edge_density(aes(fill = value)) +
    geom_edge_link(aes(width = abs(value),color=type,edge_alpha = abs(value)), alpha = 0.2,
                   arrow = arrow(length = unit(5, 'mm')), end_cap = circle(0, 'mm')) + 
    geom_node_point(size = 8,colour=rgb(220/255,20/255,60/255)) +
    geom_node_text(aes(label = name), size = 4, repel = TRUE) +
    geom_edge_loop()+
    scale_color_brewer(palette = "Set2") +
    theme_graph()
}
#plot_grn(data)


pvalues <- read.delim("E:/DRCTdb/ignore/LDSC_results/sample1/pvalues.tsv")

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
#plot_disease_heatmap(pvalues)
