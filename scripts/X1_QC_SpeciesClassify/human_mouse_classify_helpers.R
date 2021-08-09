# The purpose of this script is to contain helper functions
# for classifying the spots into data vs mouse.

library(stringr)
source('scripts/utils/helpers.R')

getSpeciesScore <- function(spatial, species_prefix='hg38', 
                            assay='Spatial', matrix='counts', sum_method=T){
  # Sums across the species genes for each spot
  # Uses gene prefixes to determine which species the genes belongs to.
  count_matrix <- slot(spatial@assays[[assay]], matrix)
  species_genes <- getSpeciesGenes(count_matrix, species_prefix)
  print('Example species genes:')
  print(species_genes[1:3])
  species_count_matrix <- count_matrix[species_genes,]
  if (sum_method) { species_scores <- apply(species_count_matrix, 2, sum) }
  else {  species_scores <- apply(species_count_matrix>0, 2, sum) }
  return( species_scores )
}

classifyByCuttoffs <- function(scores, human_cutoffs, mouse_cutoffs){
  # Sets hard-cutoffs for what's considered data versus mouse
  # based on the scores. 
  # Assumes data score in column 1 of scores, and mouse in column 1
  species <- character(nrow(scores))
  
  human_bool_x <- scores[,1]>human_cutoffs[1]
  human_bool_y <- scores[,2]<human_cutoffs[2]
  species[human_bool_x & human_bool_y] <- 'human'
  
  mouse_bool_x <- scores[,1]<mouse_cutoffs[1]
  mouse_bool_y <- scores[,2]>mouse_cutoffs[2]
  species[mouse_bool_x & mouse_bool_y] <- 'mouse'
  
  species[species==""] <- 'mix'
  
  return(species)
}

