# (29/06/2021) -> Performs clusterProfiler 
# GSEA on mSigDB using the human DE results.
# Tutorial here:
#   http://yulab-smu.top/clusterProfiler-book/chapter12.html
# 
# Complicated analysis object prevents inputting enrichment
# from other results. 
# 
# INPUT: data/DE_out/Pseudo_Limma_human/
#        de_results_PseudoLimma_human.xlsx
# 
# OUTPUT: data/DE_out/Pseudo_Limma_human/gsea_out/*
#         figure_components/DE_figures/

library(xlsx)
library(stringr)
work_dir <- rstudioapi::getSourceEditorContext()$path
work_dir <- str_split(work_dir, 'scripts/')[[1]][1]
setwd(work_dir)
source('scripts/X2_DEAnalysis/clusterProfiler_helpers.R')

data_dir <- 'data/DE_out/Pseudo_Limma_human/'
out_dir <- 'data/DE_out/Pseudo_Limma_human/gsea_out/'
out_plots <- 'figure_components/DE_figures/'

################################################################################
                 # Performing GSEA enrichment analysis #
################################################################################
# Using mSigDB as the database, see here:
#  http://yulab-smu.top/clusterProfiler-book/chapter3.html#msigdb-analysis

# 7.2.0 for msigdbr to work with vissE
library(stringr)
library(clusterProfiler)
library(enrichplot)
library(xlsx)
library(ggplot2)
library(vissE)
library(GSEABase)
library(msigdb)

species <- "human"
suffix <- "human"

######### Loading the reference gene sets ###############
msigdbr_show_species()

db_t2g <- load_db(species)

############ Loading the query gene sets ####################
file_name <- paste(data_dir, 'de_results_PseudoLimma_', suffix, '.xlsx', sep='')
geneList <- load_geneList(file_name)

############ Performing GSEA ####################
em <- GSEA(geneList, TERM2GENE = db_t2g, by='fgsea',
           minGSSize = 20, pvalueCutoff = .01, nPerm=10000)
# Adding in the percentage #
counts <- apply(em@result, 1, function(vals) {length(str_split(vals['core_enrichment'], '/')[[1]])})
perc <- round((counts / em@result[,'setSize'])*100, 1)
em@result['count'] <- counts
em@result['Percentage'] <- perc

result <- em@result
result <- result[order(result[,'NES'], decreasing = T),]
sets <- rownames(result)
print(sets[str_detect(sets, 'HALLMARK')])
print(length(sets))

##### Writing out the results #######
save_gseaResults(em, result, out_dir, suffix)

################################################################################
                      # Visualising the results #
################################################################################
file_name <- paste(out_dir, 'gsea_results_', suffix, '.rds', sep='')
em <- readRDS(file_name)
###### NOTE I flipped the NES score here so upregulation would be red to be consistent
###### in the figure, & changed the figure legend in the illustrator file to reflect this.
em@result[,'nes'] <- -em@result[,'NES']
result <- em@result
result <- result[order(result[,'NES'], decreasing = T),]

######### ClusterProfiler approach ##########
selected_terms <- c("GOBP_NEUROGENESIS", "GOBP_NEURON_DIFFERENTIATION", 
                    "GOBP_NEURON_DEVELOPMENT", "GOBP_CELL_CELL_SIGNALING",
                    "GOBP_AXON_DEVELOPMENT", "GOBP_NERVOUS_SYSTEM_PROCESS",
                    "GOBP_SYNAPTIC_SIGNALING", "GOBP_REGULATION_OF_SECRETION",
                    "HALLMARK_MYC_TARGETS_V1", "GOBP_CELL_DIVISION", "GOBP_DNA_REPLICATION",
                    "GOBP_CELL_CYCLE_PHASE_TRANSITION", "GOBP_CELL_CYCLE", "GOBP_HORMONE_TRANSPORT",
                    "HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS")
selected_terms <- format_go_names(selected_terms)

em_formatted <- em
rownames(em_formatted@result) <- format_go_names(rownames(em_formatted@result))
em_formatted@result$ID <- rownames(em_formatted@result)
em_formatted@result$Description <- rownames(em_formatted@result)
names(em_formatted@geneSets) <- rownames(em_formatted@result)

# the dotplot of the GSEA results #
p <- dotplot(em_formatted, x="GeneRatio", showCategory=selected_terms, size='Percentage',
             font.size=22) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              text=element_text(face='bold', size=20), panel.border = element_blank(),
              axis.line = element_line())
plot_name <- paste(out_plots, species, suffix, 'GSEAdotplot.pdf', sep='_')
pdf(plot_name, width=8)
print(p)
dev.off()

# The emap plot #
edo <- pairwise_termsim(em_formatted)
p <- emapplot(edo, showCategory=selected_terms, pie_scale=1, repel=T,
         cex_label_category=1.2, layout='nicely', color='nes',
         )+theme(text=element_text(face='bold', size=1))
p$layers[[3]]$aes_params[['fontface']] <- "bold"
plot_name <- paste(out_plots, species, suffix, 'GSEAemap.pdf', sep='_')
pdf(plot_name)
print(p)
dev.off()




