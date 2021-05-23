# The purpose of this script is to perform spot-specific enrichment analysis
# in order to understand the data.

# Following the instructions here on how to load in the data, except using the
# SCTransform normalised data:
# http://spatialgiotto.rc.fas.harvard.edu/giotto.visium.brain.html

################################################################################
              # Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/"
save_data_dir <- paste(work_dir, 'data/DE_out/Pseudo_TMM_Limma_Voom/giotto_out/', sep='')
save_dir <- paste(work_dir, 'figure_components/giotto_enrich_PseudoLimma_out/', sep='')
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

## Creating the Giotto object from Seurat (function implemented in helpers.R) ##
giotto <- genGiottoFromSeurat(merged,
                              save_plot=T, show_plot=T, save_dir=save_dir)

#########Loading gene sets ######
filter_genes <- giotto@gene_metadata[[1]]
load_SHH = F
load_gsea = T
load_gbm = F
min_gsea = T
gene_lists <- loadGeneSets(load_SHH = load_SHH, load_gsea = load_gsea, load_gbm=load_gbm,
                           human=T, min_gsea=min_gsea,#filter_genes=filter_genes
                           ) # defined in helpers.R
list_names <- names(gene_lists)

####### Perform the enrichment ! ########
# See details on how to run here:
# http://spatialgiotto.rc.fas.harvard.edu/giotto_guide_to_cell_type_enrichment.html

marker_matrix <- makeSignMatrixPAGE(list_names, gene_lists)

# Need to look at overlap between the sig genes and the genes we have !
enrich_results <- runPAGEEnrich(giotto, marker_matrix, expression_values='normalized',
                              p_value=T, n_times=100,
                              return_gobject = F)

# Adding the enrichment scores and p-values to the merged Seurat object #
enrich_scores <- as.data.frame( enrich_results$matrix )
rownames(enrich_scores) <- enrich_scores[,1]
enrich_scores <- enrich_scores[,!(colnames(enrich_scores) %in% c('cell_ID'))]
colnames(enrich_scores) <- str_c(colnames(enrich_scores), 'enrich scores', sep=' ')

# Replace inf values with max not inf #
for (coli in 1:ncol(enrich_scores)) {
  isinf <- is.infinite(enrich_scores[,coli])
  enrich_scores[isinf,coli] <- max(enrich_scores[isinf==F,coli])
}

keep_cols <- setdiff(colnames(merged@meta.data), colnames(enrich_scores))
merged@meta.data <- cbind(merged@meta.data[,keep_cols], enrich_scores)

# Below not really in use #
enrich_pvals <- list()
enrich_binary <- list() # binary labels, enriched or not
pval_labels <- c()
binary_labels <- c()
for (i in 1:length(list_names)) {
  list_name <- list_names[i]
  name_bool <- enrich_results$DT[['cell_type']] == list_name
  subset <- enrich_results$DT[name_bool,]
  pvals <- subset[['pval']]
  names(pvals) <- subset[['cell_ID']] 
  # Ensuring save with correct order 
  label <- paste(list_name, '.pval', sep='')
  pval_labels <- c(pval_labels, label)
  enrich_pvals[[i]] <- pvals[rownames(merged@meta.data)]
  
  binary_label <- character(length(pvals))
  names(binary_label) <- subset[['cell_ID']] 
  binary_label[pvals<.05] <- paste('enriched ', list_name, sep='')
  binary_label[pvals>=.05] <- paste('not enriched ', list_name, sep='')
  enrich_binary[[i]] <- binary_label
  binary_labels <- c(binary_labels, paste(list_name, '.binary', sep=''))
}
enrich_pvals <- as.data.frame(enrich_pvals)
colnames(enrich_pvals) <- pval_labels
enrich_binary <- as.data.frame(enrich_binary)
colnames(enrich_binary) <- binary_labels

keep_cols <- setdiff(colnames(merged@meta.data), colnames(enrich_pvals))
merged@meta.data <- cbind(merged@meta.data[,keep_cols], enrich_pvals)
keep_cols <- setdiff(colnames(merged@meta.data), colnames(enrich_binary))
merged@meta.data <- cbind(merged@meta.data[,keep_cols], enrich_binary)


# Saving the output #
out_name <- save_data_dir
run_bool <- c(load_SHH, load_gsea, load_gbm, min_gsea)
run_state <- c('SHH', 'gsea', 'gbm', 'mingsea')
for (i in 1:length(run_state)) {
  if (i==1) {middle=''} else {middle='_'}
  if (run_bool[i]) {
    out_name <- paste(out_name, middle, run_state[i], sep='')
  }
}
out_name <- paste(out_name, '_scores.txt', sep='')
write.table(enrich_scores, out_name, sep='\t', quote=F)


# Write plots to plot all of this on a consistent scale !!!
#     & using Tuan's prefferred visualisation parameters. 
print(colnames(enrich_scores))
SpatialPlot(merged, features ="SHH.C enrich scores", images='slice1')#_C1_untreated')
# SHH1 & SHH2 make no sense, but SHH.A, SHH.B, SHH.C very interesting !
# NOTE, the pvalue is awful!

### Based on the enrichment plots, binarise the enrich scores to get cell annotations #####

# Adding in the spot names from the species classification #
merged_orig <- readRDS('data/seurat_rds/merged.rds')

species <- merged_orig$species[names(merged$orig.ident)]

merged <- AddMetaData(merged, species, 'species')

# Saving the .rds !!!
saveRDS(merged, 'data/seurat_rds/all.rds')



### Old code from experimenting with how Giotto takes input !!!! #######

# raw_paths <- getRawSpatialPaths(suffix='spatial/tissue_positions_list.csv')
# 
# # Loading the spot info #
# x <- data.table::fread(raw_paths[1]) # returns a list type !
# colnames(x) = c('barcode', 'in_tissue', 'array_row', 'array_col', 'col_pxl', 'row_pxl')
# giotto_spatial = x[,.(row_pxl,-col_pxl)] # spatial locations need to simply consist 
#                                         # of row pixels and -values for the col pixels. 




