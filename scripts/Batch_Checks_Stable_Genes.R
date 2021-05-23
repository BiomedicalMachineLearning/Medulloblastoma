# The purpose of this script is to perform PCA on the sample spots to make
# sure there is not any batch effect in the data.

################################################################################
                          # Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/"
setwd(work_dir)
source('scripts/helpers.R')

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scMerge)
  library(scater)
  library("readxl")
  library(VennDiagram)
  library(stringr)
})

# Stable genes #
data(segList, package='scMerge')

human_stable <- segList$human$human_scSEG
human_stable <- str_c(rep('hg38-', length(human_stable)), human_stable)

mouse_stable <- segList$mouse$mouse_scSEG
mouse_stable <- str_c(rep('mm10-', length(mouse_stable)), mouse_stable)

stable <- c(human_stable, mouse_stable)

################################################################################
                      # Loading the data #
################################################################################
merged <- readRDS('data/seurat_rds/all.rds')

mouse_stable <- intersect(mouse_stable, getSpeciesGenes(merged@assays$SCT@counts, 'mm10-'))
human_stable <- intersect(human_stable, getSpeciesGenes(merged@assays$SCT@counts, 'hg38-'))
stable <- c(human_stable, mouse_stable)

Idents(merged) <- 'species'
species <- unique(merged$species)

for (i in 1:length(species)) {
  spec_i <- species[i]
  spec <- subset(merged, idents=c(species[i]))
  if (spec_i=='human') {stable_i <- human_stable}
  else if (spec_i=='mouse') {stable_i <- mouse_stable}
  else {stable_i <- stable}
  
  dims <- c(RunPCA, RunUMAP)
  dim_names <- c('pca', 'umap')
  for (j in 1:2){
    dim_i <- dims[[j]]
    dim_name <- dim_names[j]
    
    # Using the stable genes #
    spec <- dim_i(spec, feature=stable_i)
    plot_title <- paste(spec_i, 'spots')
    plot <- DimPlot(spec, reduction=dim_name, group.by = 'sample') + ggtitle(plot_title)
    plot_file <- paste('figure_components/batch_check/stable_genes_', dim_name, spec_i, '.pdf', sep='')
    dealWithPlot(T, plot_file, plot)
    
    # Using hvgs #
    spec <- dim_i(spec, features=merged@assays$SCT@var.features)
    plot_title <- paste(spec_i, 'spots')
    plot <- DimPlot(spec, reduction=dim_name, group.by = 'sample') + ggtitle(plot_title)
    plot_file <- paste('figure_components/batch_check/hvgs_', dim_name, spec_i, '.pdf', sep='')
    dealWithPlot(T, plot_file, plot)
  }
}








