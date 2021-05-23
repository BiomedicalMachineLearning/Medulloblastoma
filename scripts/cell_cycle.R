# Purpose of this script is to score the cells on cell cycle #

################################################################################
                            # Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/"
setwd(work_dir)
source('scripts/helpers.R')

library(Seurat)
library(ggplot2)
library(stringr)

################################################################################
                  # Loading in the data #
################################################################################
all <- readRDS('data/seurat_rds/all.rds')

DefaultAssay(all) <- "SCT"

################################################################################
                  # Cell cycle Scoring #
################################################################################
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
s.genes <- str_c(rep("hg38-", length(s.genes)), s.genes)
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- str_c(rep("hg38-", length(g2m.genes)), g2m.genes)

all <- CellCycleScoring(all, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(all[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(all, features = c("hg38-PCNA", "hg38-TOP2A", "hg38-MCM6", "hg38-MKI67"), ncol = 2)

# Visualising the cell cycle #
all <- RunPCA(all, features = c(s.genes, g2m.genes))
DimHeatmap(all, dims = c(1, 2))

all <- RunPCA(all)
all <- RunUMAP(all, reduction = "pca", dims = 1:30, n.components=2L)
DimPlot(all)

DimPlot(all, reduction = "umap", #split.by='treatment', 
        group.by = 'Phase', pt.size=1, dims=c(1,2))

FeaturePlot(all, reduction = "umap", #split.by='treatment', 
        features = 'hg38-ATOH1', pt.size=1, dims=c(1,2))



