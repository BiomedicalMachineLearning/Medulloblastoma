################################################################################
                      # Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/"
setwd(work_dir)
source('scripts/helpers.R')

library(Seurat)
library(ggplot2)

################################################################################
                    # Loading the data #
################################################################################
path_info <- getSpatialPaths()
data_paths <- path_info[[1]]
data_names <- path_info[[2]]

spatials <- readRDSFiles(data_paths, data_names)

################################################################################
            # Plotting gene counts for each species #
################################################################################
source('scripts/human_mouse_classify_helpers.R')

out_dir <- 'figure_components/species_plots/'

spatial_species <- list()
for (i in 1:length(spatials)){
  # Calculating the species scores #
  human_genes <- getSpeciesScore(spatials[[i]], species_prefix='hg38', 
                                  assay='SCT', matrix='data', sum_method=F) 
  mouse_genes <- getSpeciesScore(spatials[[i]], species_prefix='mm10',
                                  assay='SCT', matrix='data', sum_method=F) 
  
  scores = as.data.frame(cbind(human_genes, mouse_genes))
  rownames(scores) <- names(spatials[[i]]$orig.ident)
  colnames(scores) <- c('human_genes', 'mouse_genes')
  
  spatials[[i]] <- AddMetaData(spatials[[i]], scores)
  
  features <- c('human_genes', 'mouse_genes')
  for (j in 1:2) {
    feature <- features[j]
    score_plot <- SpatialFeaturePlot(spatials[[i]], 
                                     features=feature, pt.size.factor=1)
    plot_name <- paste(out_dir, data_names[i], '_', feature, '.pdf', sep='')
    dealWithPlot(save_plot=T, plot_name, score_plot, width=8)
  }
}


################################################################################
              # Predicting human vs. mouse tissue spots #
################################################################################
source('scripts/human_mouse_classify_helpers.R')

out_dir <- 'figure_components/species_plots/'
save_plot <- T

human_cuts <- list(c(250, 400), c(250, 400), c(250, 300), c(250, 250))
mouse_cuts <- list(c(250, 250), c(250, 250), c(250, 250), c(250, 250))

spatial_species <- list()
for (i in 1:length(spatials)){
  # Calculating the species scores #
  human_scores <- getSpeciesScore(spatials[[i]], species_prefix='hg38', 
                                  assay='SCT', matrix='data', sum_method=F) 
  mouse_scores <- getSpeciesScore(spatials[[i]], species_prefix='mm10',
                                  assay='SCT', matrix='data', sum_method=F) 
  
  scores = as.data.frame(cbind(human_scores, mouse_scores))
  rownames(scores) <- names(spatials[[i]]$orig.ident)
  
  human_cutoffs <- human_cuts[[i]]
  mouse_cutoffs <- mouse_cuts[[i]]
  
  species <- classifyByCuttoffs(scores, human_cutoffs, mouse_cutoffs)
  scores[,'species'] <- species
  spatial_species[[names(spatials)[i]]] <- scores
  
  spatials[[i]] <- AddMetaData(spatials[[i]], scores)
  
  plot <- ggplot(scores, aes(x=mouse_scores, y=human_scores, color=species)) + geom_point()
  plot_name <- paste(out_dir, data_names[i], '_speciesCounts.pdf', sep='')
  dealWithPlot(save_plot, plot_name, plot)
  
  normalPlot('species', spatials[[i]], data_names[i], save_plot=save_plot)
  
  # Saving the meta data #
  species_meta <- as.data.frame(spatials[[i]]$species)
  colnames(species_meta) <- c('species')
  write.table(species_meta, paste('data/spot_meta/', names(spatials)[i], 
                                  '_speciesCounts.txt', sep=''), quote=F)
}

##### Saving the results ########
saveRDSFiles(data_paths, spatials)







