# Purpose of this script is to create the GSEA figures;
# namely, to visualise the statistics for each terms & the overlap of the genes of each term,
# in both a network & a dotplot. 

# Following the example provided here:
# http://yulab-smu.top/clusterProfiler-book/chapter12.html 

# Good paper on benchmarking enrichment methods here:
# https://academic.oup.com/bib/article/22/1/545/5722384

# OUTPUT: figure_components/PseudoLimma/DE_figures/*

################################################################################
                    # Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/"
setwd(work_dir)
source('scripts/helpers.R')

library(ggplot2)
library(stringr)
library(hexbin)
library(cowplot)
library(stringr)

# Dependencies for GSEA result visualisations #
library(DOSE)
library(enrichplot)

# Running through the example data provided to get a better idea of structure #
data(geneList)
de <- names(geneList)[abs(geneList) > 2]
edo <- enrichDGN(de) # Main datastructure on which the other DE functions #
colnames(edo@result)
# "ID"          "Description" "GeneRatio"   "BgRatio"     "pvalue"      "p.adjust"    "qvalue"      "geneID"      "Count"

library(reticulate)

pd <- import("pandas")

out_dir <- 'figure_components/PseudoLimma/DE_figures/'

################################################################################
                    # Reading in the data #
################################################################################
gsea <- pd$read_excel(
          'data/DE_out/Pseudo_TMM_Limma_Voom/gsea_out/interesting_results/gsea_summary_human.treatment.xlsx',
                      engine='openpyxl')

termToGenes <- list()
for (i in 1:length(gsea$Term)) {
  termToGenes[[gsea$Term[i]]] <- str_split(gsea$genes[1], ';')[[1]]
}

# TODO 






