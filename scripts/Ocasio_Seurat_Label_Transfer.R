# The purpose of this script is to use the 
# labels stored from the SHH inhibitor paper to label our 
# spatial Medulloblastoma data.


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

# Third party data #
ocasio <- readRDS('data/third_party_data/Ocasio2019_NatCom_scRNA/GT.rds')
misc <- as.data.frame(ocasio@misc)
ocasio@misc <- NULL
ocasio <- UpdateSeuratObject(ocasio)











