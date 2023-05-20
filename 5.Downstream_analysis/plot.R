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




