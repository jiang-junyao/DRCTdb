run_scenic_main <- function(cis_dbDir = "E:\\public\\cistarget\\human cistarget ranking",
                            nCores=1){
  library(SCENIC)
  scenicOptions <- initializeScenic(org="hgnc", 
                                    dbDir="E:\\public\\cistarget\\human cistarget ranking", nCores=1)
  rna_test = as.matrix(rna_test)
  genesKept <- geneFiltering(rna_test, scenicOptions)
  exprMat_filtered <- rna_test[genesKept, ]
  runCorrelation(exprMat_filtered, scenicOptions)
  exprMat_filtered_log <- log2(exprMat_filtered+1) 
  runGenie3(exprMat_filtered_log, scenicOptions)
  exprMat_log <- log2(rna_test+1)
  scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget"))
}