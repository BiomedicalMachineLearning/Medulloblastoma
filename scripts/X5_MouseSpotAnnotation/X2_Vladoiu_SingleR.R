# Performs the automated annotation of the mouse spots using the Vladoiu mouse
#                                                 cerebellum data as a reference.
#
#    INPUT: * data/scanpy_h5ads/*_all_species_SME.h5ad
#           * data/third_party_data/
#                      Vladiou2019_Nature_scRNA/
#                                Vladoiu2019_lateRef.h5ad
#
#    OUTPUT: * data/spot_meta/species_classify_v2/
#                       *Vladoiu_singleR_scores_mouse.txt
#            * figure_components/MouseAnnot_figures/

######################################################################
                    # Environment Setup #
######################################################################
library(stringr)

work_dir <- rstudioapi::getSourceEditorContext()$path
work_dir <- str_split(work_dir, 'scripts/')[[1]][1]
setwd(work_dir)

source('scripts/utils/helpers.R')
source('scripts/utils/SingleR_helpers.R')

# Set to your conda environment if you're using one #
#Sys.setenv(RETICULATE_PYTHON = "~/opt/miniconda3/envs/MB/bin/python/")
library(SingleR)
library(Seurat)
library(ggplot2)
library(reticulate)
ad <- import('anndata')
source_python('scripts/utils/st_helpers_min.py')

data_dir <- 'data/third_party_data/Vladiou2019_Nature_scRNA/'
data_dir2 <- 'data/scanpy_h5ads/'
out_dir <- 'data/spot_meta/species_classify_v2/'
out_plots <- 'figure_components/MouseAnnot_figures/'

samples <- c('A1', 'B1', 'C1', 'D1')
seur_names <- c('A1_treated_', 'B1_treated_', 'C1_untreated_', 'D1_untreated_')

######################################################################
                  # Loading in the datasets #
######################################################################
mouse_mix <- list()
for (i in 1:length(samples)) {
  ad_ <- ad$read_h5ad(paste0(data_dir2, samples[i], '_all_species_SME.h5ad'))
  mouse_ad <- species_split(ad_, species='mouse')
  mouse_mix[[samples[i]]] <- t(mouse_ad$to_df())
}

vlad_ad <- ad$read_h5ad(paste0(data_dir, 'Vladoiu2019_lateRef.h5ad'))
vlad_expr <- t(vlad_ad$to_df())
vlad_labels <- as.character(vlad_ad$obs$cell_labels)

common_genes <- rownames(vlad_expr)
for (i in 1:length(mouse_mix)) {
  common_genes <- intersect(common_genes, rownames(mouse_mix[[i]]))
}

vlad_expr_mouse <- vlad_expr[common_genes,]

######################################################################
           # Running through the SingleR pipeline for mouse #
######################################################################
vlad_aggr <- aggregateReference(vlad_expr_mouse, vlad_labels, ncenters=3)
trained <- trainSingleR(vlad_aggr, vlad_aggr$label, 
                        de.n=200, de.method='t')
print('Genes for correlation:')
print(length(trained$common.genes))

# For Mouse #
singleR_map(mouse_mix, trained,
            '_mouse', out_plots, out_dir, prefix='Vladoiu')






