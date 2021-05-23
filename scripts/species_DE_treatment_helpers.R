# The purpose of this script is to impliment helper functions for 
# calling the DE for the spatial data !

getOrderedResults <- function(de_results, padj_cutoff){
  # Using a DESeq2 object, pulls out the DE results. Extra arguments parsed to results()#
  de_results <- de_results[order(de_results$padj),]
  de_results <- de_results[which(de_results$padj<padj_cutoff),]
  return (de_results)
}

directionSeperate <- function(de_results) {
  # Separates into up vs down regulated genes #
  upReg = de_results[which(de_results$log2FoldChange>0),]
  downReg = de_results[which(de_results$log2FoldChange<0),]
  
  newFormat = list()
  newFormat[["up"]] <- upReg
  newFormat[["down"]] <- downReg
  return(newFormat)
}

formatResults <- function(de_results, padj_cutoff, ...){
  # Gets the DESeq2 results in a desired format #
  de_results <- getOrderedResults(de_results, padj_cutoff)
  de_results <- directionSeperate(de_results)
  return(de_results)
}


