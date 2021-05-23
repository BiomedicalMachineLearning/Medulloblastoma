# Purpose of this script is to load in the factor loadings from the cell2loc output
# to visualise how these look spatially.

################################################################################
                      # Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/"
setwd(work_dir)
source('scripts/helpers.R')

library(Seurat)
library(ggplot2)
library(stringr)

################################################################################
                        # Loading the data #
################################################################################
merged <- readRDS('data/seurat_rds/all.rds')

factor_scores <- read.table('data/decon_out/labelTransfer_human_vladoiu_cell2LocScores.txt',
                              sep='\t', row.names=1, header=T)

fact_names <- str_c(rep("labelTransfer_factors", ncol(factor_scores)),
                    as.character(1:ncol(factor_scores)))
colnames(factor_scores) <- fact_names

merged <- AddMetaData(merged, factor_scores)

SpatialPlot(merged, features='labelTransfer_factors10', images = 'slice1')

# NOTES: * fact1 interesting
#        * fact2 noise
#        * fact3 mostly in mouse..
#        * fact4 intersting
#        * fact8 in mixed samples..
#        * fact8 very interesting !
# Almost all factors only relevant to mouse spots, need to redo removing these !!!!
# also need to redo the label transfer while removing these. 

