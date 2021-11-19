# Using cluster profiler to do an over-representation analysis on
#                the mouse LR pairs with atleast 10 significant spots on the
#                tumour border. But this time adding them all to the same gene
#                list, instead of performing independently for each sample.
#
#         INPUT:  * data/cci/mouse/interface_overlaps.xlsx
#         OUTPUT: * figure_components/cci_figures/_mouse_LR-genes_GSEAemap.pdf
#                 * data/cci/mouse/gsea_results_mouse_LR-genes.xlsx

################################################################################
                          # Environment setup #
################################################################################
Sys.setenv(RETICULATE_PYTHON = "~/opt/miniconda3/envs/MB/bin/python/")
library(reticulate)
sc <- import('scanpy')
pd <- import('pandas')

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

data_dir <- 'data/cci/mouse/'
out_dir <- data_dir
out_plots <- 'figure_components/cci_figures/'

samples <- c('A1', 'B1', 'C1', 'D1')
species <- 'mouse'
suffix <- 'LR-genes'

################################################################################
                  # Loading & formatting data #
################################################################################
### Loading the ranked LRs consistent across samples ###
overlap_df <- pd$read_excel(paste0(data_dir, 'interface_overlaps.xlsx'),
                            index_col=as.integer(0))
final_lrs <- overlap_df[,'A1_B1_C1_D1']

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
em <- enrichGO(genes, org.Mm.eg.db, ont='BP', keyType='SYMBOL', 
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











