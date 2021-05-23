# The purpose of this script is to explore getting 
# Seurat spatial up and running with the Medulla Blastoma data.

# Going through the Seurat Spatial Vignette located here:
# https://satijalab.org/seurat/articles/spatial_vignette.html 

################################################################################
                      # Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/"
data_dir <- paste(work_dir, 'data/Visium8_A1_Hybrid_treated', sep='')
setwd(work_dir)

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

################################################################################
                          # Loading in the data #
################################################################################
spatial <- Load10X_Spatial(data_dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
)

# Some basic QC plotting #
plot1 <- VlnPlot(spatial, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(spatial, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

########################## SCTransform normalisation #############################
# NOTE: in the vignette, it includes a comparison between this and log-transform #
spatial <- SCTransform(spatial, assay = "Spatial", verbose = FALSE)

# Plotting gene expression in the Seurat space #
SpatialFeaturePlot(spatial, features = c("hg38-HPCA", "hg38-OAZ3"))
# Great !!!!
# TODO:
#   * save the data as a .rds file
#   * add in command-line input so that will automatically generate a publication
#     worthy pdf of inputted gene names for each slice. 







