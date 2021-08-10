# Purpose of this script is to use house functions which are generally
# useful across scripts.

library(Seurat)
library(stringr)
#library(rhdf5)

work_dir <- rstudioapi::getSourceEditorContext()$path
work_dir <- str_split(work_dir, 'scripts/')[[1]][1]

################################################################################
                      # Data loading functions #
################################################################################
getRawSpatialPaths <- function(prefix='', suffix=''){
  # Gets the paths to the raw data from 10X #
  raw_paths <- c('Visium8_A1_Hybrid_treated/', 'Visium8_B1_Hybrid_treated/', 
                 'Visium8_C1_Hybrid_untreated/', 'Visium8_D1_Hybrid_untreated/')
  data_paths <- c()
  for (i in 1:length(raw_paths)){
    data_paths <- c(data_paths, paste(prefix, 'data/', raw_paths[i], suffix, sep=''))
  }
  return(data_paths)
}

getSpatialPaths <- function(){
  # Gets the default paths to the data #
  data_names <- c('A1_treated', 'B1_treated', 'C1_untreated', 'D1_untreated')
  data_paths <- c()
  for (i in 1:length(data_names)){
    data_paths <- c(data_paths, paste('data/seurat_rds/', 
                                      data_names[i], '.rds', sep=''))
  }
  
  out <- list('data_paths'=data_paths, 'data_names'=data_names)
  return(out)
}

readRDSFiles <- function(data_paths, data_names) {
  spatials <- list()
  for (i in 1:length(data_paths)){
    spatials[[data_names[i]]] <- readRDS(data_paths[i])
  }
  return(spatials)
}

saveRDSFiles <- function(data_paths, data_objects){
  for (i in 1:length(data_paths)){
    saveRDS(data_objects[[i]], file=data_paths[i])
  }
}

replaceByMapping <- function(mapping, labels){
  # Mapping is a list with the values in labels as keys, and the 
  # desired label to replace the labels with as values.
  new_labels <- character(length(labels))
  for (i in 1:length(mapping)){
    new_labels[labels==names(mapping)[i]] <- mapping[[i]]
  }
  
  return(new_labels)
}

loadGeneSets <- function(load_gsea_v2=T, load_SHH=T,
                         human=T, filter_genes=NULL, add_prefix=F) {
  # Loads the compiled interesting gene sets based on my own 
  # GSEA analysis and also from third-party papers
  # optional arguments indicate which gene sets to load.

  ######## Loading the gsea_v2 results using clusterProfiler ###########
  if (load_gsea_v2) {
    file_name <- "data/DE_out/Pseudo_Limma_human/gsea_out/gsea_results_human.xlsx"
    gsea_humantreat <- as.data.frame(read_excel(file_name,
                                                sheet = "human_human"))
    rownames(gsea_humantreat) <- gsea_humantreat[,'ID']
    terms_v2 <- c("HALLMARK_E2F_TARGETS")
    gsea_humantreat <- gsea_humantreat[terms_v2,]

    terms_list_v2 <- list()
    for (i in 1:length(terms_v2)) {
      genes <- str_split(gsea_humantreat[i, 'core_enrichment'], '/')[[1]]
      terms_list_v2[[i]] <- genes
    }
  } else {terms_v2 <- c(); terms_list_v2 <- list()}

  ####### Load in the reference genes for different SHH subtypes ######
  sign_names <- c()
  sign_list <- list()
  if (load_SHH) {
    sign_names <- c()
    sign_list <- list()
    file_dir <- 'data/third_party_data/'
    file_names <- c("Hovestadt_SHHA-cellcycle_program.txt",
                    "Hovestadt_SHHB-undiff_program.txt",
                    "Hovestadt_SHHC-neuron_program.txt")
    for (i in 1:length(file_names)) {
      df <- read.table(paste(file_dir, file_names[i], sep=''), header=T)
      sign_names <- c(sign_names, colnames(df))
      sign_list[[i]] <- df[,1]
    }
  }

  # Compiling the gene lists to run #
  list_names <- c(sign_names, terms_v2)
  gene_lists <- c(sign_list, terms_list_v2)
  names(gene_lists) <- list_names
  
  # Making sure is in the desired format & filtered to detected genes #
  for (i in 1:length(gene_lists)) {
    list_name <- list_names[i]
    # Adding prefix # 
    if (human & add_prefix) {prefix <- 'hg38-'} 
    else if (!human & add_prefix) {prefix <- 'mm10-'}
    else {prefix <- ''}
    genes <- str_c(rep(prefix, length(gene_lists[[list_name]])), gene_lists[[list_name]])
    
    final_genes <- c()
    for (j in 1:length(genes)) {
      gene <- genes[j]
      if (is.null(filter_genes) || #No filter genes, keep everything.
         (!is.null(filter_genes) && gene %in% filter_genes)) {
        final_genes <- c(final_genes, gene)
      }
    }
    gene_lists[[list_names[i]]] <- final_genes
  }
  
  return(gene_lists)
}

loadhdf5data <- function(h5File) {
  # Function for loading data that was saved from pandas as a R dataframe #
  # NOTE, ended up using reticulate instead to load using pandas, more flexible #
  
  listing <- h5ls(h5File)
  # Find all data nodes, values are stored in *_values and corresponding column
  # titles in *_items
  data_nodes <- grep("_values", listing$name)
  name_nodes <- grep("_items", listing$name)
  
  data_paths = paste(listing$group[data_nodes], listing$name[data_nodes], sep = "/")
  name_paths = paste(listing$group[name_nodes], listing$name[name_nodes], sep = "/")
  
  columns = list()
  for (idx in seq(data_paths)) {
    data <- data.frame(t(h5read(h5File, data_paths[idx])))
    names <- t(h5read(h5File, name_paths[idx]))
    entry <- data.frame(data)
    colnames(entry) <- names
    columns <- append(columns, entry)
  }
  
  data <- data.frame(columns)
  
  return(data)
}

################################################################################
                  # Data manipulation functions #
################################################################################
mergeData <- function(spatials, data_names){
  # Merges the Seurat data #
  
  spatial1 <- spatials[[1]]
  other_spatials <- list()
  for (i in 2:length(spatials)) {
    other_spatials[[i-1]] <- spatials[[i]]
  }
  merged <- merge(spatial1, other_spatials, add.cell.ids=data_names)
  
  return(merged)
}

genGiottoFromSeurat <- function(merged, human_genes=T,
                                save_plot=T, show_plot=T, save_dir='.', 
                                python_path=NULL, assay='SCT', ...){
  # Takes in a Seurat object and converts to a Giotto object #
  # -> Assumes Seurat object SCT normalised 
  # -> Assumes only desired data genes, and hence only subset to these
  # -> Extra arguments parsed to 'createGiottoObject' function.
  print("Reminder: this function subset to just the human genes !")
  raw_expr <- merged@assays[[assay]]@counts
  norm_expr <- merged@assays[[assay]]@data
  #scale_expr <- merged@assays$SCT@scale.data
  if (human_genes) {
    prefix <- 'hg38-'
  } else {
    prefix <- 'mm10-' 
  }
  
  genes <- getSpeciesGenes(raw_expr, prefix)
  raw_expr <- raw_expr[genes,]
  norm_expr <- norm_expr[genes,]
  
  gene_names <- str_replace(genes, prefix, "")
  gene_names <- toupper(gene_names)
  rownames(raw_expr) <- gene_names
  rownames(norm_expr) <- gene_names
  
  cell_meta <- merged@meta.data
  
  # Creating the spatial information to keep giotto happy (but dosn't actually make sense)
  image_info <- merged@images
  coords <- image_info[[1]]@coordinates
  for (i in 2:length(image_info)) {
    coords <- rbind(coords, image_info[[i]]@coordinates)
  }
  # Making sure is in-line #
  print(all(rownames(coords) == colnames(norm_expr)))
  
  # Re-arrange the coords to satisfy Giotto input, as indicated in the above url
  giotto_coords <- coords[,c('imagerow', 'imagecol')]
  colnames(giotto_coords) <- c('row_pxl', 'col_pxl')
  giotto_coords[,2] <- -giotto_coords[,2]
  
  # Create giotto object, NOTE this is using brad-site conda environment for python
  # How to create giotto object details here:
  # https://rdrr.io/github/RubD/Giotto/man/createGiottoObject.html
  myinst=createGiottoInstructions(save_plot=T, show_plot=T, 
                                  save_dir = save_dir, python_path=python_path)
  giotto <- createGiottoObject(raw_exprs = raw_expr, 
                               norm_expr = norm_expr,
                               #norm_scaled_expr = scale_expr,
                               spatial_locs = giotto_coords, 
                               instructions = myinst, 
                               cell_metadata = cell_meta)
  
  return(giotto)
}


################################################################################
                    # Gene query functions #
################################################################################

getSpeciesGenes <- function(count_matrix, species_prefix, strip_prefix=F){
  # Gets the genes with the corresponding species_prefix #
  gene_names <- rownames(count_matrix)
  species_bool <- str_detect(gene_names, species_prefix)
  species_genes <- gene_names[species_bool]
  
  if (strip_prefix) {
    species_genes <- str_replace(species_genes, species_prefix, "")
  }
  
  return(species_genes)
}

getSpeciesEquivalentGenes <- function(genes, count_matrix, 
                                      prefixes=c('hg38-', 'mm10-')){
  # Given a list of genes specified in normal data format (MITF),
  # retrieve the genes which match the inputted genes, assuming
  # species genes specified with prefixes: 'hg38-', 'mm10-'
  gene_by_species <- list()
  for (i in 1:length(prefixes)) {
    species_genes <- getSpeciesGenes(count_matrix, prefixes[i])
    split <- str_split(species_genes, prefixes[i])
    
    match_genes <- c()
    for (j in 1:length(split)) {
      formatted <- toupper(split[[j]][2])
      if (formatted %in% genes){
        match_genes <- c(match_genes, species_genes[j])
      }
    }
    
    gene_by_species[[prefixes[i]]] <- match_genes
  }
  
  return( gene_by_species )
}

convertSpeciesGenes <- function(gene_names, prefix='hg38-') {
   # Converts to the other species genes from the genes given
   # using the current species indicated in the prefix.
   genes <- str_replace(gene_names, prefix, "")
   new_genes <- c()
   for (i in 1:length(genes)) {
      if (prefix=='hg38-') {
         first <- substring(genes[i], 1, 1)
         rest <- tolower(substring(genes[i], 2))
         new_gene <- paste(first, rest, sep='')
      
      } else {
         new_gene <- toupper(genes[i])    
      }
      new_genes <- c(new_genes, new_gene)
   } 

   return( new_genes )
}

calc_cell_RLE <-
  function (expr_mat, spikes = NULL) 
    # Calculating the RLE values #
  {
    RLE_gene <- function(x) {
      #if (median(unlist(x)) > 0) {
        log((x + 1)/(median(unlist(x)) + 1))/log(2)
      #}
      #else {
      #  rep(NA, times = length(x))
      #}
    }
    if (!is.null(spikes)) {
      RLE_matrix <- t(apply(expr_mat[-spikes, ], 1, RLE_gene))
    }
    else {
      RLE_matrix <- t(apply(expr_mat, 1, RLE_gene))
    }
    #cell_RLE <- apply(RLE_matrix, 2, median, na.rm = T)
    return(RLE_matrix)
  }

################################################################################
                        # Plotting functions #
################################################################################
dealWithPlot <- function(save_plot, plot_name, plot, ...) {
  # Deals with the plot !
  if (save_plot) {
    print('Saved plot to:')
    print(plot_name)
    ggsave(plot_name, plot, ...)
  } else { print(plot)  }
}

spatialExprPlot <- function(genes, spatial, data_name, save_plot = T,
                            out_dir = 'figure_components/expr_plots/',
                            pt.size=1, ...) {
  # Plots the spatial expression for a given seurat object #
  count_matrix <- spatial@assays$SCT@counts
  genes_by_species <- getSpeciesEquivalentGenes(genes, count_matrix)
  
  # Plotting the genes for each section and saving output #
  for (i in 1:length(genes_by_species)) {
    species_genes <- genes_by_species[[i]]
    for (j in 1:length(species_genes)){
      gene <- species_genes[j]
      plot <- SpatialFeaturePlot(spatial, features = c(gene), 
                                 pt.size.factor = pt.size, ...)
      plot_name <- paste(out_dir, data_name, '-', 
                         species_genes[j], '.pdf', sep='')
      
      # Dealing with output #
      dealWithPlot(save_plot, plot_name, plot)
    }
  }
}

normalPlot <- function(info, spatial, data_name, save_plot=T,
                       out_dir = 'figure_components/', 
                       pt.size=1, ...) {
  # Performing normal plotting #
  meta_info <- names(spatial)
  gene_info <- names(spatial@assays$SCT@counts)
  # if (info %in% meta_info) {
  plot <- SpatialPlot(spatial, group.by=info, 
                        pt.size.factor=pt.size, ...)
  # } else{
  #   plot <- SpatialPlot(spatial, group.by=info, 
  #                       pt.size.factor=pt.siz, ...)
  # }
  
  plot_name <- paste(out_dir, as.character(info), 
                     '-', data_name, '.pdf', sep='')
  dealWithPlot(save_plot, plot_name, plot)
}














