# Loads in the fetal brain data, & uses this as a
# reference to label cells by dominant cell type.
# First download the fetal human brain data to match INPUT:
# Ref. paper: https://www.nature.com/articles/s41586-020-2157-4
# Download link: https://db.cngb.org/HCL/gallery.html?tissue=Fetal-Brain3
#
#    INPUT: * data/scanpy_h5ads/*_all_species_SME.h5ad
#           * data/third_party_data/HCL2020/Fetal-Brain3_dge.txt
#           * data/third_party_data/HCL2020/Fetal-Brain3_Anno.txt
#
#    OUTPUT: * data/spot_meta/*FetalBrain3singleR_scores.txt
#            * figure_components/HumanAnnot_figures/

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

data_dir <- 'data/third_party_data/HCL2020/'
data_dir2 <- 'data/scanpy_h5ads/'
out_dir <- 'data/spot_meta/species_classify_v2/'
out_plots <- 'figure_components/HumanAnnot_figures/'

samples <- c('A1', 'B1', 'C1', 'D1')
seur_names <- c('A1_treated_', 'B1_treated_', 'C1_untreated_', 'D1_untreated_')

######################################################################
                  # Loading in the datasets #
######################################################################
#### Preprocessing the data! #####
# Loading in the reference data & creating Seurat object #
ref_counts <- read.csv(paste0(data_dir, 'Fetal-Brain3_dge.txt'), row.names=1)
ref_meta <- read.csv(paste0(data_dir, 'Fetal-Brain3_Anno.csv'), row.names=1)
print(all(rownames(ref_meta)==colnames(ref_counts)))
ref_seurat <- CreateSeuratObject(ref_counts, meta.data=ref_meta)

# Normalising with sctransform & getting normalised data #
ref_seurat <- SCTransform(ref_seurat, return.only.var.genes = F)
# Performing UMAP #
ref_seurat <- RunPCA(ref_seurat, features = VariableFeatures(object = ref_seurat))
ref_seurat <- RunUMAP(ref_seurat, dims=1:50)

DimPlot(ref_seurat, group.by='CT')

## Saving the UMAP coordinates for visualisation in scanpy #
ref_umap <- as.data.frame(ref_seurat@reductions$umap@cell.embeddings)
write.table(ref_umap, paste0(data_dir, 'Fetal-Brain3_seuratUMAP.txt'),
            sep='\t')

## Pulling out the normalised data to run with SingleR #
ref_expr <- as.matrix(ref_seurat@assays$SCT@data)
ref_labels <- as.character(ref_seurat$CT)

###### Loading the query Visium data #####
human_mix <- list()
for (i in 1:length(samples)) {
  ad_ <- ad$read_h5ad(paste0(data_dir2, samples[i], '_all_species_SME.h5ad'))
  human_ad <- species_split(ad_, species='human')
  human_mix[[samples[i]]] <- t(human_ad$to_df())
}

common_genes <- rownames(ref_expr)
for (i in 1:length(human_mix)) {
  common_genes <- intersect(common_genes, rownames(human_mix[[i]]))
}
print(length(common_genes)) # 10678

ref_expr_human <- ref_expr[common_genes,]

######################################################################
      # Running through the SingleR pipeline for human #
######################################################################
ref_aggr_h <- aggregateReference(ref_expr_human, ref_labels, ncenters=3)
trained_h <- trainSingleR(ref_aggr_h, ref_aggr_h$label, 
                        de.n=350, de.method='classic')
print('Genes for correlation:')
print(length(trained_h$common.genes)) #3338

# For Human #
singleR_map(human_mix, trained_h, '_human',
            out_plots, out_dir, prefix='FetalBrain3')





