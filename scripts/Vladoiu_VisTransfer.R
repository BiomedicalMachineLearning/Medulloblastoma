# The purpose of this script is to visualise the results from
# the label transfer procedure. 

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
merged <- readRDS('data/seurat_rds/all.rds')

################################################################################
                # Results from the Vladoiu transfer #
################################################################################
#### Functions related to analysis below #####
labelByMaxProb <- function(seurat_obj, assay.name='human_vladoiu_preds') {
  # Generates cell labels based on max probs and adds as extra cell meta data #
  cell_probs <- seurat_obj@assays[[assay.name]]@data
  max_probs <- cell_probs['max',]
  spot_labels <- c()
  for (i in 1:length(max_probs)) {
    spot_probs <- cell_probs[,i]
    cell <- names(spot_probs)[which(spot_probs==max_probs[i])[1]]
    spot_labels <- c(spot_labels, paste(assay.name, cell, sep="_"))
  }
  return(spot_labels)
}

###################### HUMAN ###############################
####### Adding labels for the data genes ######
DefaultAssay(merged) <- "human_vladoiu_preds"
human_pred_labels <- labelByMaxProb(merged)
merged <- AddMetaData(merged, human_pred_labels, col.name='human_vladoiu_pred_labels')

pred_labels <- unique(human_pred_labels)
human_pred_freqs <- list()
for (i in 1:length(pred_labels)) {
  human_pred_freqs[[pred_labels[i]]] <- length(which(human_pred_labels==pred_labels[i]))
}

####### Using data/mix spots ########
Idents(merged) <- 'species'
human <- subset(merged, idents=c('human', 'mix'))

Idents(human) <- 'human_vladoiu_pred_labels'
SpatialPlot(human, group.by='human_vladoiu_pred_labels', images='slice1_B1_treated',
            pt.size.factor=1, interactive=T)

###################### MOUSE ###############################
####### Adding labels for the mouse genes ######
DefaultAssay(merged) <- "mouse_vladoiu_preds"
mouse_pred_labels <- labelByMaxProb(merged, 'mouse_vladoiu_preds')
merged <- AddMetaData(merged, mouse_pred_labels, col.name='mouse_vladoiu_pred_labels')

pred_labels <- unique(mouse_pred_labels)
mouse_pred_freqs <- list()
for (i in 1:length(pred_labels)) {
  mouse_pred_freqs[[pred_labels[i]]] <- length(which(mouse_pred_labels==pred_labels[i]))
}

####### Using mouse spots ########
Idents(merged) <- 'species'
mouse <- subset(merged, idents=c('mouse'))

Idents(mouse) <- 'mouse_vladoiu_pred_labels'
SpatialPlot(mouse, group.by='mouse_vladoiu_pred_labels', images='slice1',
            pt.size.factor=1, interactive=T)


# TODO cluster based on the probabilities... use NMF approach of cell2location which
#       takes as input the proportion dataframe;
#       https://cell2location.readthedocs.io/en/latest/cell2location.downstream_models.html#module-cell2location.models.CoLocatedGroupsSklearnNMF






