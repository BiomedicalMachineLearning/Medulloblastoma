# The purpose of this script is to load in the 
# spatial seurat data, and then generate plots of the spatial 
# expression and save to figure_components/expr_plots

################################################################################
                    # Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/"
setwd(work_dir)
source('scripts/helpers.R')

library(Seurat)
library(ggplot2)

################################################################################
                    # Loading the data #
################################################################################
path_info <- getSpatialPaths()
data_paths <- path_info[[1]]
data_names <- path_info[[2]]

spatials <- readRDSFiles(data_paths, data_names)

################################################################################
                  # Generating the spatial plots #
################################################################################
# Getting the equivalent genes from both species #
# genes <- c(#'CD44', 'SERPINE1', 'TNFRSF1A',
#            'CDK1', 'CCNB2', 'CKS1B', 'CDKN3', 'NASP', # Cell cycle markers/E2F targets
#            'FOXM1', 'VIM',  'PTTG1', #SP1 targets
#            'NR2F1', 'LRP8', #VEGF A UP.V1 DN
#            'FABP3', 'SPARC', 'CRYAB', 'CKB', 'NCAM1', 'MYL6B', 'ITGB1', 'MYH7', # Myogenesis
#            'CALM1', 'CALM2', 'TGFB3', 'TGFBR1', 'VDAC1', 'VDAC2' # Cell Senescence
#            )
genes <- c('IL4', 'IGF1', 'IGF1R', 'IGFBP3') # some CCI expression
for (i in 1:length(spatials)){
  spatialExprPlot(genes, spatials[[i]], data_names[i], save_plot = T,
                  pt.size=1, alpha=.7)
}


################################################################################
              # Generating the violin plots #
################################################################################
merged <- readRDS('data/seurat_rds/merged.rds')
Idents(merged) <- 'species'
human <- subset(merged, idents=c('human'))

Idents(human) <- 'treatment'

for (i in 1:length(genes)) {
  gene_name <- paste('hg38-', genes[i], sep='')
  plot <- VlnPlot(object = human, features = gene_name, split.by = 'treatment')
  plot_name <- paste('figure_components/violin_plots/', gene_name, 
                     '_humanspots_treatvcontrol.pdf', sep='')
  print(plot_name)
  dealWithPlot(T, plot_name, plot)
}

################################################################################
                # Generating the dot plots #
################################################################################
genes_wPrefix <- str_c('hg38-', genes)
plot <- DotPlot(human, features = genes_wPrefix) + RotatedAxis()
plot_name <- 'figure_components/dot_plots/pathway_markers.pdf'
dealWithPlot(T, plot_name, plot, width=15, height=4)


# gene_groups <- list('Cell cycle markers/E2F targets'=c('CDK1', 'CCNB2', 'CKS1B', 'CDKN3', 'NASP'))
# for (i in 1:length(gene_groups)) {
#   group_genes <- gene_groups[[i]]
#   genes_wPrefix <- str_c('hg38-', genes)
#   DotPlot(data, features = genes_wPrefix) + RotatedAxis()
# }






