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
merged <- readRDS('data/seurat_rds/all.rds')

out_dir <- 'figure_components/PseudoLimma/'
giotto_out <- paste(out_dir, 'giotto_enrich_plots/', sep='')
violin_out <- paste(out_dir, 'violin_plots/', sep='')

dir.create(out_dir, showWarnings = FALSE)
dir.create(giotto_out, showWarnings = FALSE)
dir.create(violin_out, showWarnings = FALSE)

################################################################################
              # Generating the enrichment plots #
################################################################################
Idents(merged) <- 'species'
human <- subset(merged, idents=c('human', 'mix'))

enrich_scores_bool <- str_detect(colnames(merged@meta.data), 'enrich scores')
enrich_score_names <- colnames(merged@meta.data)[enrich_scores_bool]
#enrich_score_names <- c('gbm_bulk enrich scores', 'gbm_emt enrich scores')
image_names <- names(merged@images)
image_to_exp <- list('slice1'= 'A1_treated', 'slice1_B1_treated' = 'B1_treated',
                      'slice1_C1_untreated' = 'C1_untreated', 
                      'slice1_D1_untreated'= 'D1_untreated')
#pt_sizes <- c(2, 1.5, 1, 1)
pt_sizes <- c(2.5, 1.5, 1.5, 2)
for (i in 1:length(enrich_score_names)) {
  enrich_score_name <- enrich_score_names[i]
  max_cutoff <- max( merged[[enrich_score_name]] )
  for (j in 1:length(image_names)) {
    image <- image_names[j]
    plot <- SpatialPlot(human, features=enrich_score_name, images=image,
                        pt.size=pt_sizes[j], alpha=.8,
                        max.cutoff = max_cutoff)
    plot <- plot + scale_fill_gradientn(#limits = c(min_cutoff,max_cutoff),
                                       colours=c("blue4", "yellow"))
    plot_name <- paste(giotto_out, image_to_exp[[image]], '_', 
                       enrich_score_name, '.pdf', sep='')
    dealWithPlot(T, plot_name, plot)
  }
}

# Plotting also with violin plots #
# Idents(merged) <- 'species'
# data <- subset(merged, idents=c('data', 'mix'))

Idents(human) <- 'treatment'

for (i in 1:length(enrich_score_names)) {
  enrich_score_name <- enrich_score_names[i]
  plot <- VlnPlot(object = human, features = enrich_score_name, 
                  split.by = 'treatment')
  plot_name <- paste(giotto_out, 'violin_', 
                     enrich_score_name,
                     '_humanspots_treatvcontrol.pdf', sep='')
  print(plot_name)
  dealWithPlot(T, plot_name, plot)
}

####### NOTE below not used yet !!!, these are from copying the other vis script
################################################################################
                  # Generating the spatial plots #
################################################################################
# Getting the equivalent genes from both species #
genes <- c(#'CD44', 'SERPINE1', 'TNFRSF1A',
           'CDK1', 'CCNB2', 'CKS1B', 'CDKN3', 'NASP', # Cell cycle markers/E2F targets
           'FOXM1', 'VIM',  'PTTG1', #SP1 targets
           'NR2F1', 'LRP8', #VEGF A UP.V1 DN
           'FABP3', 'SPARC', 'CRYAB', 'CKB', 'NCAM1', 'MYL6B', 'ITGB1', 'MYH7', # Myogenesis
           'CALM1', 'CALM2', 'TGFB3', 'TGFBR1', 'VDAC1', 'VDAC2' # Cell Senescence
           )
for (i in 1:length(spatials)){
  spatialExprPlot(genes, spatials[[i]], data_names[i], save_plot = F,
                  pt.size=1, alpha=.7)
}


################################################################################
              # Generating the violin plots #
################################################################################
Idents(merged) <- 'species'
human <- subset(merged, idents=c('human'))

Idents(human) <- 'treatment'

for (i in 1:length(genes)) {
  gene_name <- paste('hg38-', genes[i], sep='')
  plot <- VlnPlot(object = human, features = gene_name, split.by = 'treatment')
  plot_name <- paste(violin_out, gene_name, 
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






