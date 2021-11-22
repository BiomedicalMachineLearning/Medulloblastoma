# Using cluster profiler to do an over-representation analysis on
#                the LR pairs with atleast 4 significant spots on the tumour
#                border.
#
#                INPUT:  * data/cci/human/*_border_enriched_LRs.txt
#                OUTPUT: * figure_components/cci_figures/
#                                                   _human_LR-genes_GSEAemap.pdf
#                        * data/cci/mouse/gsea_results_human_LR-genes.xlsx

################################################################################
                          # Environment setup #
################################################################################
Sys.setenv(RETICULATE_PYTHON = "~/opt/miniconda3/envs/MB/bin/python/")
library(reticulate)
sc <- import('scanpy')

library(xlsx)
library(stringr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(vissE)
library(GSEABase)
library(msigdb)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(DOSE)
source('scripts/X6_CellCellInteraction/clusterProfiler_helpers.R')

work_dir <- rstudioapi::getSourceEditorContext()$path
work_dir <- str_split(work_dir, 'scripts/')[[1]][1]
setwd(work_dir)

data_dir <- 'data/cci/human/'
out_dir <- data_dir
out_plots <- 'figure_components/cci_figures/'

samples <- c('A1', 'B1', 'C1', 'D1')
species <- 'human'
suffix <- 'LR-genes'

################################################################################
                  # Loading & formatting data #
################################################################################
### Loading the ranked LR lists ###
lr_dfs <- list()
for (i in 1:length(samples)){
  samp <- samples[i]
  path <- paste0(data_dir, samp, '_border_enriched_LRs.txt')
  lr_df <- read.table(path)
  lr_dfs[[samp]] <- lr_df
}

### Getting all lrs used across samples ###
all_lrs <- c()
for (i in 1:length(samples)) {
  lr_df <- lr_dfs[[i]]
  plot(1:dim(lr_df)[1], lr_df[,"border_spot_counts"])
  lrs <- rownames(lr_df)
  lrs <- str_remove_all(lrs, 'hg38-')
  all_lrs <- c(all_lrs, lrs)
  #lrs_split <- as.data.frame( str_split(lrs, '_') )
  #ls <- as.character(lrs_split[1,])
  #rs <- as.character(lrs_split[2,])
  #genes <- unique( c(ls, rs) )
  #print(length(genes))
  #samp_genes[[samples[i]]] <- genes
}

### Getting those present in two samples ####
all_lrs <- unique( all_lrs )
lr_counts <- integer(length(all_lrs))
for (i in 1:length(samples)){
  lr_df <- lr_dfs[[i]]
  lrs <- rownames(lr_df)
  lrs <- str_remove_all(lrs, 'hg38-')
  for (j in 1:length(lrs)){
    lr_counts[all_lrs==lrs[j]] <- lr_counts[all_lrs==lrs[j]]+1
  }
}

#### Retrieving the genes ###
final_lrs <- all_lrs[lr_counts>1]

lrs_split <- as.data.frame( str_split(final_lrs, '_') )
ls <- as.character(lrs_split[1,])
rs <- as.character(lrs_split[2,])
genes <- unique( c(ls, rs) )

################################################################################
                          # Performing ORA #
################################################################################
######### Loading the reference gene sets ###############
msigdbr_show_species()

db_t2g <- load_db(species)

################ Performing the ORA ##########################
em <- enrichGO(genes, org.Hs.eg.db, ont='BP', keyType='SYMBOL', 
                 pvalueCutoff=.01, qvalueCutoff=.05)

result <- em@result
sig_results <- result[,'p.adjust']<.01
result <- result[sig_results,]

##### Writing out the results #######
save_gseaResults(em, result, out_dir, paste(species, suffix, sep='_'))

#### Visualising the results ######
edo <- pairwise_termsim(em)
p <- emapplot(edo, showCategory=28, 
              pie_scale=1, repel=T,
              cex_label_category=1.2, layout='nicely', color='p.adjust',
)+theme(text=element_text(face='bold', size=1))
p$layers[[3]]$aes_params[['fontface']] <- "bold"
plot_name <- paste(out_plots, species, suffix, 'GSEAemap.pdf', sep='_')
pdf(plot_name, width=12, height=10)
print(p)
dev.off()











