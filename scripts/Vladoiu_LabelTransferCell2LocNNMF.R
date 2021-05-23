# Created: 24th of March 2021
# Edit:

# The purpose of this script is to perform  

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



