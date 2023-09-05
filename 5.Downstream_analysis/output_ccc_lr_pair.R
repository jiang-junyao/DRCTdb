library(CellChat)
setwd('E:\\DRCTdb\\ignore\\downstream_result')
dir1 = dir()
for (i in dir1) {
  ccc = readRDS(paste0(i,'/cellchat_ccc.rds'))
  for (j in 1:length(ccc)) {
    if (!is.null(ccc[[j]])) {
      if (!is.na(ccc[[j]])) {
        ct_ccc = ccc[[j]]
        if (sum(ct_ccc@net[["count"]])>0) {
          groupSize = as.numeric(table(ct_ccc@idents))
          # filenames2 = paste0(i,'/ccc/',names(ccc)[j],'.svg')
          # svg(filename = filenames2, width = 4, height = 4)
          # print(netVisual_circle(ct_ccc@net$weight, vertex.weight = groupSize,
          #                        weight.scale = T, label.edge= F,
          #                        title.name = "",vertex.label.cex = 0.7))
          # dev.off()
          netVisual_circle(ct_ccc@net$weight, vertex.weight = groupSize,
                           weight.scale = T, label.edge= F,
                           title.name = "",vertex.label.cex = 0.7)
          ct_lr = netVisual_bubble(ct_ccc)[['data']]
          write.table(ct_lr,paste0(i,'/ccc/',names(ccc)[j],'.txt'),quote = F,
                      row.names = F,sep = '\t')
        }

      }

    }

  }
}
