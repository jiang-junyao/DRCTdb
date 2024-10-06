run_cellchat <- function(a1,metadata,group,species,
                         search = c('Cell-Cell Contact','ECM-Receptor',
                                    'Secreted Signaling')){
  library(CellChat)
  cellchat <- createCellChat(object = a1, meta = metadata, group.by = group)
  if (species =='mm') {
    db = CellChatDB.mouse
  }
  if (species =='hs') {
    db = CellChatDB.human
  }
  if (species =='zf') {
    db = CellChatDB.zebrafish
  }
  CellChatDB <- db 
  CellChatDB.use <- subsetDB(CellChatDB, search = search)
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  future::plan("multicore", workers = 4)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  return(cellchat)
}
