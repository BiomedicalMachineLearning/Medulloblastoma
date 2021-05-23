# Created: 22nd March 2021
# Edit: 

# The purpose of this script is to create a Seurat v4 object using the Vladiou 
# mouse medulloblastoma scRNA-seq data, in preparation for SCTransform 
# normalisation and subsequent label transfer onto the spatial data.

# INPUT: data/third_party_data/Vladoiu2019_Nature_scRNA/
#             Vladoiu2019_cell_counts.hdf, Vladoiu2019_cell_meta.txt
# OUTPUT: data/third_party_data/Vladoiu2019_Nature_scRNA/
#             vladoiu.rds

################################################################################
                        # Environment setup #
################################################################################
#data_dir <- 'data/third_party_data/Vladoiu2019_Nature_scRNA/'
data_dir <- '/30days/uqbbalde/MedullaBlastoma/Vladoiu2019_Nature_scRNA/'

# Pointing to the conda python for reticulate #
Sys.setenv(RETICULATE_PYTHON = "/90days/uqbbalde/.conda/envs/MedullaBlastoma/bin/python3")

library(Seurat)
library(ggplot2)
library(reticulate)
pandas <- import('pandas')

################################################################################
                          # Loading the data #
################################################################################
data_path <- paste(data_dir, 'Vladoiu2019_cell_counts.hdf', sep='')
counts <- pandas$read_hdf(data_path)

meta_path <- paste(data_dir, 'Vladoiu2019_cell_meta.txt', sep='')
cell_meta <- read.table(meta_path, header=T, sep='\t', row.names=1)

#### Creating the seurat object !!! #####
seurat_ <- CreateSeuratObject(counts, meta.data=cell_meta)

################################################################################
               # Performing SCTransform normalisation #
################################################################################
seurat_ <- SCTransform(seurat_, assay = "RNA", verbose = T,
            		return.only.var.genes=F #Ensures all genes returned !
			)

# Dim reduction #
seurat_ <- RunPCA(seurat_)
seurat_ <- RunUMAP(seurat_, dims=1:30)

# Saving the rds output #
saveRDS(seurat_, paste(data_dir, 'vladoiu2019.rds', sep=''))














