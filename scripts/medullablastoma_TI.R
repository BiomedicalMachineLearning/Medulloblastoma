# The purpose of this script is to attempt trajectory inference on the data to 
# see how this corresponds to the previous patterns identified via the GSEA
# analysis. 

################################################################################
            # Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/"
setwd(work_dir)
source('scripts/helpers.R')

library(Seurat)
library(ggplot2)
library(stringr)
library(slingshot)

################################################################################
                    # Loading the data #
################################################################################
path_info <- getSpatialPaths()
data_paths <- path_info[[1]]
data_names <- path_info[[2]]

spatials <- readRDSFiles(data_paths, data_names)

# Determining set of genes which are consistently variable across experiments #
var_genes <- VariableFeatures(spatials[[1]])
for (i in 2:length(spatials)) {
  var_genes <- setdiff(var_genes, VariableFeatures(spatials[[i]]))
}
var_genes_bool <- str_detect(var_genes, 'hg38-')
var_genes <- var_genes[var_genes_bool]

merged <- readRDS('data/seurat_rds/merged.rds')

## Subsetting to just the data/mixed spots ##
Idents(merged) <- 'species'
human <- subset(merged, idents=c('human', 'mix'))
VariableFeatures(human) <- var_genes

# TODO try clustering based on the enrichment scores!!! #

human <- RunPCA(human, assay = "SCT", verbose = FALSE)
human <- FindNeighbors(human, reduction = "pca", dims = 1:30)
human <- FindClusters(human, verbose = FALSE)
human <- RunUMAP(human, reduction = "pca", dims = 1:30)

DimPlot(human, reduction = "umap", label = TRUE, group.by='sample')
SpatialPlot(human, images = 'slice1_C1_untreated', 
            features='SHH.C enrich scores')




