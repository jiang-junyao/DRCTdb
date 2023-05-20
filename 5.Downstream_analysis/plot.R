data = read.table('../ignore/downstream_result/sample1/grn_cor04/vCM_Atrial_fibrillation.txt',header = T)
data$type = ifelse(data$value>0,'positive','negative')
plot_grn <- function(data){
  library(ggraph)
  ggraph(data,layout = 'fr') + 
    geom_edge_density(aes(fill = value)) +
    geom_edge_link(aes(width = abs(value),color=type,edge_alpha = abs(value)), alpha = 0.2) + 
    geom_node_point(size = 8,colour=rgb(220/255,20/255,60/255)) +
    geom_node_text(aes(label = name), size = 4, repel = TRUE) +
    scale_color_brewer(palette = "Set2") +
    theme_graph()
}
plot_grn(data)


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
plot_disease_heatmap(pvalues)
