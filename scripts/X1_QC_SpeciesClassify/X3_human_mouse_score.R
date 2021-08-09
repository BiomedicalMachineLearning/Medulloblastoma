#Sums to sctransformed normalised gene expression for human/mouse genes to
# generate human/mouse scores per spot in each sample, to subsequently allow
# for species clasification.
#
#                    INPUT: * seurat_rds/*treated.rds
#                    OUTPUT: * data/spot_meta/*_species.txt

################################################################################
                      # Environment setup #
################################################################################
library(Seurat)
library(ggplot2)
library(stringr)

work_dir <- rstudioapi::getSourceEditorContext()$path
work_dir <- str_split(work_dir, 'scripts/')[[1]][1]
setwd(work_dir)
source('scripts/utils/helpers.R')
source('scripts/X1_QC_SpeciesClassify/human_mouse_classify_helpers.R')

out_dir <- 'data/spot_meta/'

################################################################################
                      # Loading the data #
################################################################################
path_info <- getSpatialPaths()
data_paths <- path_info[[1]]
data_names <- path_info[[2]]

spatials <- readRDSFiles(data_paths, data_names)

################################################################################
                  # Getting human/mouse scores #
################################################################################
spatial_species <- list()
for (i in 1:length(spatials)){
  # Calculating the species scores #
  human_scores <- getSpeciesScore(spatials[[i]], species_prefix='hg38', 
                                  assay='SCT', matrix='data') 
  mouse_scores <- getSpeciesScore(spatials[[i]], species_prefix='mm10',
                                  assay='SCT', matrix='data') 
  
  scores = as.data.frame(cbind(human_scores, mouse_scores))
  rownames(scores) <- names(spatials[[i]]$orig.ident)

  # Saving the meta data #
  write.table(scores, 
              paste0('data/spot_meta/', names(spatials)[i], '_species.txt'), 
              quote=F)
}




