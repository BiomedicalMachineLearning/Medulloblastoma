# The purpose of this script is to load in each of the tissue slices,
# normalise using SCTransform (separately), annotated into data vs mouse spots,
# and then save to rds files.
# NOTE: this includes Ryan's filtering !!!
# INPUT: ryan_ids/*, VISIUM*/
# OUTPUT: seurat_rds/*.rds

################################################################################
                              # Environment setup #
################################################################################
work_dir <- "/30days/uqbbalde/MedullaBlastoma/"
setwd(work_dir)

data_dir <- '/QRISdata/Q1851/Quan/Visium/Visium8/'

#### Constructing the paths to the spatial data #####
data_suffixes <- c('Visium8_A1_Hybrid/',
                   'Visium8_B1_Hybrid/',
                   'Visium8_C1_Hybrid/',
                   'Visium8_D1_Hybrid/')
data_names <- c('A1_treated', 'B1_treated', 'C1_untreated', 'D1_untreated')
data_dirs <- c()
for (i in 1:length(data_suffixes)) {
  data_diri <- paste(data_dir, data_suffixes[i], sep='')
  data_dirs <- c(data_dirs, data_diri)
}

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)

################################################################################
                # Using Ryan's pre-filtering of the data #
################################################################################
# Loading in Ryan's IDs after filtering #
id_dir <- 'data/ryan_ids/'
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
  print("This many cells filtered by Ryan: ")
  print(length(ids[[i]]))
  
  spatial <- subset(spatial, cells=ids[[i]])
  
  print("Length after subsetting to Ryan Ids:")
  print(length(spatial$orig.ident))
  
  # SCTransforming #
  spatial <- SCTransform(spatial, assay = "Spatial", verbose = FALSE,
                         return.only.var.genes=F #Ensures all genes returned !
                         )
  
  spatials[[data_names[i]]] <- spatial
  
  saveRDS(spatial, paste('data/seurat_rds/', data_names[i], '.rds', sep=''))
}

# Merging and saving the data #
merged <- mergeData(spatials, data_names)
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

saveRDS(merged, 'data/seurat_rds/merged.rds')

################################################################################
          # Using the data by themselves, applying SCTransform #
################################################################################
# Loading in Ryan's IDs after filtering #
id_dir <- 'data/ryan_ids/'
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





