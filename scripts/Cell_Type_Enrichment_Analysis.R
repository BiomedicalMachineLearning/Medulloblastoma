# Purpose of this script is to perform enrichment 
# for cell type specific genes per spot in an attempt to label
# spots in terms of their active pathways. 
# INPUT: data/seurat_rds/integrated.rds

################################################################################
                      # Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/"
setwd(work_dir)
source('scripts/helpers.R')

library(Seurat)
library(ggplot2)
library(stringr)
library(AUCell)

################################################################################
                      # Loading in the data #
################################################################################
human.combined <- readRDS('data/seurat_rds/integrated.rds')

exprMatrix <- as.matrix(human.combined@assays$Spatial@counts)

###### Filtering out bad genes ###### 
# Ribosome/MT genes #
all_genes <- rownames(exprMatrix)
annoying_genes <- str_detect(all_genes, pattern = '-RPL')
annoying_genes <- annoying_genes | str_detect(all_genes, pattern = '-RPS')
annoying_genes <- annoying_genes | str_detect(all_genes, pattern = '-MT-') 
annoying_genes <- annoying_genes | str_detect(all_genes, pattern = '-Rpl')
annoying_genes <- annoying_genes | str_detect(all_genes, pattern = '-Rps')
annoying_genes <- annoying_genes | str_detect(all_genes, pattern = '-mt-') 
remaining_genes <- all_genes[annoying_genes==F]

human_genes_bool <- str_detect(remaining_genes, pattern='hg38-')
human_genes <- remaining_genes[human_genes_bool]

exprMatrix <- exprMatrix[human_genes,]

# Getting genes detected... #
gene_counts <- apply(exprMatrix>0, 2, sum)

more_genes <- gene_counts>1000
  
human.combined <- AddMetaData(human.combined, gene_counts, 'gene_counts')
human.combined <- AddMetaData(human.combined, more_genes, 'high_gene_spots')

FeaturePlot(human.combined, "gene_counts",
            pt.size = 1, order=T)

DimPlot(human.combined, group.by='high_gene_spots', pt.size = 1, shuffle=T)

exprMatrix <- exprMatrix[,more_genes]

################################################################################
                # Now performing the AUC Scoring #
################################################################################
cells_rankings <- AUCell_buildRankings(exprMatrix)

###### Loading in marker genes for the different Vladoiu datasets ###### 
cell_type_genes <- read.csv('data/third_party_data/Vladoiu2019_Nature_scRNA/vladoiu_de_genes.txt',
                            sep='\t', row.names=1)
gene_sets <- as.list(cell_type_genes)
for (i in 1:length(gene_sets)){
  n_len <- 400
  gene_sets[[i]] <- str_c( rep("hg38-", n_len), toupper(gene_sets[[i]][1:n_len]) )
}

cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings, aucMaxRank=200)

hist(cells_AUC@assays@data$AUC[2,])

i <- 20
human.combined <- AddMetaData(human.combined, 
                              cells_AUC@assays@data$AUC[i,], names(gene_sets)[i])

FeaturePlot(human.combined, names(gene_sets)[i], pt.size = 1, order=T)




















