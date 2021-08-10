# clusterProfiler_helpers.R (29/06/2021) -> Helper function for 
#                                           clusterProfiler_[species].R.

library(msigdbr)
library(xlsx)
library(stringr)
library(dplyr)

load_db <- function(species){
  # Loads the reference gene sets #
  species_ <- list('human'= 'Homo sapiens', 'mouse'= 'Mus musculus')
  species_prefix <- list('human'= 'hg38-', 'mouse'= 'mm10-')
  
  cats <- c('C5', 'H')
  subcats <- c('GO:BP', NA)
  db <- msigdbr(species = species_[[species]], category=cats[1], subcategory=subcats[1])
  for (i in 2:length(cats)) {
    if (!is.na(subcats[i])) {subcat <- subcats[i]} else {subcat<-NULL}
    db <- full_join(db, msigdbr(species = species_[[species]], 
                                category=cats[i], subcategory=subcat))
  }
  
  db_t2g <- db %>% dplyr::select(gs_name, gene_symbol)
  return(db_t2g)
}

load_geneList <- function(file_name, prefix='hg38-') {
  # Loads the DE gene list #
  up_df <- read.xlsx(file_name, 1)
  down_df <- read.xlsx(file_name, 2)
  
  geneList <- c(up_df[,'t.val'], down_df[,'t.val'])
  gene_names <- c(up_df[,1], down_df[,1])
  species_bool <- str_detect(gene_names, prefix)
  geneList <- geneList[species_bool]
  gene_names <- gene_names[species_bool]
  gene_names <- str_sub(gene_names, 6)
  names(geneList) <- gene_names
  geneList <- geneList[order(geneList, decreasing=T)]
  return(geneList)
}

save_gseaResults <- function(em, result, out_dir, suffix) {
  # Save the gsea results #
  file_name <- paste(out_dir, 'gsea_results_', suffix, '.rds', sep='')
  saveRDS(em, file_name)
  
  file_name <- paste(out_dir, 'gsea_results_', suffix, '.xlsx', sep='')
  sheet_name <- paste(suffix, species, sep='_')
  write.xlsx(result, file_name,
             sheetName=sheet_name, col.names=TRUE, row.names=TRUE, append=F)
}

format_go_names <- function(go_names) {
  go_names <- str_replace_all(go_names, 'GOBP_', '')
  go_names <- str_replace_all(go_names, 'HALLMARK_', '')
  go_names <- str_replace_all(go_names, '_', ' ')
  go_names <- tolower(go_names)
  return(go_names)
}

format_result <- function(em){
  result <- em@result
  result <- result[order(result[,'NES'], decreasing = T),]
  return(result)
}






