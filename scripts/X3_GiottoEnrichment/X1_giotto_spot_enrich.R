# Using the Giotto package for PAGE per-spot enrichment analysis:
#   https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02286-2#Sec2
#   http://spatialgiotto.rc.fas.harvard.edu/giotto.install.native.html#s_macos
# 
# To perform enrichment for DE E2F targets & different 
# MB-SHH transcriptional programs, which were downloaded from
# Hovestadt, et al. into data/third_party_data/.
# 
# INPUT: * data/seurat_rds/all.rds
#        * data/Pseudo_Limma_human/gsea_out/
#                         gsea_results_human.xlsx
#        * data/third_party_data/*_program.txt
# 
# OUTPUT: * data/giotto_out/*
#         * figure_components/giotto_figures/*

################################################################################
                        # Environment setup #
################################################################################
library(Seurat)
library(ggplot2)
library(Giotto)
library("readxl")
library(stringr)

work_dir <- rstudioapi::getSourceEditorContext()$path
work_dir <- str_split(work_dir, 'scripts/')[[1]][1]
setwd(work_dir)

save_data_dir <- 'data/giotto_out/'
save_dir <- 'figure_components/giotto_figures/'
setwd(work_dir)
source('scripts/utils/helpers.R')

################################################################################
                        # Loading the data #
################################################################################
merged <- readRDS('data/seurat_rds/all.rds')

## Creating the Giotto object from Seurat (function implemented in helpers.R) ##
giotto <- genGiottoFromSeurat(merged,
                              save_plot=T, show_plot=T, save_dir=save_dir)

######### Loading gene sets ######
filter_genes <- giotto@gene_metadata[[1]]
load_SHH = T
load_gsea_v2 = T
gene_lists <- loadGeneSets(load_SHH = load_SHH, load_gsea_v2 = load_gsea_v2,
                           human=T) # defined in helpers.R
list_names <- names(gene_lists)
print(gene_lists) # Looks good!

####### Perform the enrichment ! ########
# See details on how to run here:
# http://spatialgiotto.rc.fas.harvard.edu/giotto_guide_to_cell_type_enrichment.html

marker_matrix <- makeSignMatrixPAGE(list_names, gene_lists)

# Need to look at overlap between the sig genes and the genes we have !
enrich_results <- runPAGEEnrich(giotto, marker_matrix, expression_values='normalized',
                              p_value=T, n_times=100,
                              return_gobject = F)

# Retrieving the per-spot enrichment scores #
enrich_scores <- as.data.frame( enrich_results$matrix )
rownames(enrich_scores) <- enrich_scores[,1]
enrich_scores <- enrich_scores[,!(colnames(enrich_scores) %in% c('cell_ID'))]
colnames(enrich_scores) <- str_c(colnames(enrich_scores), 'enrich scores', sep=' ')

# Replace inf values with max not inf #
for (coli in 1:ncol(enrich_scores)) {
  isinf <- is.infinite(enrich_scores[,coli])
  enrich_scores[isinf,coli] <- max(enrich_scores[isinf==F,coli])
}

# Saving the output #
out_name <- save_data_dir
run_bool <- c(load_SHH, load_gsea_v2)
run_state <- c('SHH', 'gsea_v2')
for (i in 1:length(run_state)) {
  if (i==1) {middle=''} else {middle='_'}
  if (run_bool[i]) {
    out_name <- paste(out_name, middle, run_state[i], sep='')
  }
}
out_name <- paste(out_name, '_scores.txt', sep='')
write.table(enrich_scores, out_name, sep='\t', quote=F)


