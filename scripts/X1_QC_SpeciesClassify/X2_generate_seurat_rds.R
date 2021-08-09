# Generates seurat rds with sctransform normalised data which is later used
# for species scoring & giotto enrichment analysis.
#
#            INPUT: * data/Visium8_{sample_name}_Hybrid/*
#                   * data/filter_ids/{sample_name}_filtered.txt
#
#            OUTPUT: * seurat_rds/*

################################################################################
                              # Environment setup #
################################################################################
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)

work_dir <- rstudioapi::getSourceEditorContext()$path
work_dir <- str_split(work_dir, 'scripts/')[[1]][1]
setwd(work_dir)

data_dir <- paste0(work_dir, 'data/')

#### Constructing the paths to the spatial data #####
data_dirs <- c('data/Visium8_A1_Hybrid_treated/',
                   'data/Visium8_B1_Hybrid_treated/',
                   'data/Visium8_C1_Hybrid_untreated/',
                   'data/Visium8_D1_Hybrid_untreated/')
data_names <- c('A1_treated', 'B1_treated', 'C1_untreated', 'D1_untreated')

################################################################################
                # Using Ryan's pre-filtering of the data #
################################################################################
# Loading in IDs after filtering #
id_dir <- 'data/filter_ids/'
id_files <- list.files(id_dir)
id_files <- id_files[order(id_files)]
ids <- list()
for (i in 1:length(id_files)) {
  ids[[i]] <- read.table(paste(id_dir, id_files[i], sep=''))[,1]
}

spatials <- list() # Creating the spatial objects #
for (i in 1:length(data_names)) {
  # Loading the spatial data #
  spatial <- Load10X_Spatial(data_dirs[i],
                             filename = "filtered_feature_bc_matrix.h5",
                             assay = "Spatial",
                             slice = "slice1",
                             filter.matrix = TRUE,
                             to.upper = FALSE,
  )
  
  # Performing the subset to the spots Ryan is working with #
  print("This many cells loaded from seurat: ")
  print(length(spatial$orig.ident))
  print("This many cells filtered: ")
  print(length(ids[[i]]))
  
  spatial <- subset(spatial, cells=ids[[i]])
  
  print("Length after subsetting to filtered Ids:")
  print(length(spatial$orig.ident))
  
  # SCTransforming #
  spatial <- SCTransform(spatial, assay = "Spatial", verbose = FALSE,
                         return.only.var.genes=F #Ensures all genes returned !
                         )
  
  spatials[[data_names[i]]] <- spatial
  
  saveRDS(spatial, paste('data/seurat_rds/', data_names[i], '.rds', sep=''))
}

################################################################################
          # Using the data by themselves, applying SCTransform #
################################################################################
# Loading in IDs after filtering #
id_dir <- 'data/filter_ids/'
id_files <- list.files(id_dir)
id_files <- id_files[order(id_files)]
ids <- list()
for (i in 1:length(id_files)) {
  ids[[i]] <- read.table(paste(id_dir, id_files[i], sep=''))[,1]
}

spatials <- list() # Creating the spatial objects #
for (i in 1:length(data_names)) {
  # Loading the spatial data #
  spatial <- Load10X_Spatial(data_dirs[i])
  
  # Performing the subset to the spots Ryan is working with #
  print("This many cells loaded from seurat: ")
  print(length(spatial$orig.ident))
  print("This many cells filtered by Ryan: ")
  print(length(ids[[i]]))
  
  spatial <- subset(spatial, cells=ids[[i]])
  
  print("Length after subsetting to Ryan Ids:")
  print(length(spatial$orig.ident))
  
  spatials[[data_names[i]]] <- spatial
}

# Merging  #
merged <- merge(spatials[[1]], 
                y= spatials[-1], add.cell.ids=data_names)
merged <- SCTransform(merged, assay = "Spatial", verbose = T,
            return.only.var.genes=F #Ensures all genes returned !
)

# Adding meta data and saving #
name_split <- str_split(names(merged$orig.ident), '_')
sample_names <- c()
treatments <- c()
for (i in 1:length(name_split)) {
  sample_name <- paste(name_split[[i]][1:2], collapse='_')
  sample_names <- c(sample_names, sample_name)
  
  treatment <- name_split[[i]][2]
  treatments <- c(treatments, treatment)
}
merged <- AddMetaData(merged, sample_names, 'sample')
merged <- AddMetaData(merged, treatments, 'treatment')
print(unique(merged$sample))
print(unique(merged$treatment))

saveRDS(merged, 'data/seurat_rds/all.rds')





