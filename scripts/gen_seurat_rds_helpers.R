# The purpose of this script is to impliment functions that help
# generate_seurat_rds.R.

library(stringr)

getSpeciesGenes <- function(count_matrix, species_prefix){
  gene_names <- rownames(count_matrix)
  species_bool <- str_detect(gene_names, species_prefix)
  species_genes <- gene_names[species_bool]
  
  return(species_genes)
}

getSpeciesScore <- function(spatial, species_prefix='hg38', 
                            assay='Spatial', matrix='counts'){
  # Sums across the species genes for each spot
  # Uses gene prefixes to determine which species the genes belongs to.
  count_matrix <- slot(spatial@assays[[assay]], matrix)
  species_genes <- getSpeciesGenes(count_matrix, species_prefix)
  print('Example species genes:')
  print(species_genes[1:3])
  species_count_matrix <- count_matrix[species_genes,]
  species_scores <- apply(species_count_matrix, 2, sum)
  return( species_scores )
}



