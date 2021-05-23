# Purpose of this script is to simply load the data & save the enrichment scores
# to a file so can be loaded in python as well for visualisation.

################################################################################
              # Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/"
save_dir <- paste(work_dir,
                  'figure_components/giotto_enrich_PseudoLimma_out/', sep='')
setwd(work_dir)
source('scripts/helpers.R')

library(Seurat)
library(ggplot2)
library(Giotto)
library("readxl")
library(stringr)

################################################################################
                        # Loading the data #
################################################################################
merged <- readRDS('data/seurat_rds/all.rds')

score_cols <- c()
col_names <- colnames(merged@meta.data)
for (i in 1:length(col_names)) {
  if ('scores' %in% col_names[i]) {
    score_cols <- c(score_cols, col_names[i])
  }
}














