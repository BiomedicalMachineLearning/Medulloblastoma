# The purpose of this script is to derive a joint space for the 
# medulloblastoma data in order to justify the clustering performed thus far.

################################################################################
                          # Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/"
setwd(work_dir)
source('scripts/helpers.R')

library(Seurat)
library(ggplot2)
library(stringr)
library(edgeR) # Contains limma

################################################################################
                          # Loading in the data #
################################################################################
all <- readRDS('data/seurat_rds/all.rds')

# Ribosome/MT genes #
#all_genes <- rownames(all@assays$Spatial@counts)
#annoying_genes <- str_detect(all_genes, pattern = '-RPL')
#annoying_genes <- annoying_genes | str_detect(all_genes, pattern = '-RPS')
#annoying_genes <- annoying_genes | str_detect(all_genes, pattern = '-MT-') 
#annoying_genes <- annoying_genes | str_detect(all_genes, pattern = '-Rpl')
#annoying_genes <- annoying_genes | str_detect(all_genes, pattern = '-Rps')
#annoying_genes <- annoying_genes | str_detect(all_genes, pattern = '-mt-') 
#remaining_genes <- all_genes[annoying_genes==F]

#human_genes <- remaining_genes[str_detect(remaining_genes, pattern='hg38-')]

# Removing ANNOYING GENES!!! #
#all <- subset(all, features=human_genes)

# Subseting to huma/mix spots #
Idents(all) <- 'species'
human <- subset(all, idents=c('human', 'mix'))

# Adding in the cluster labels #
cluster_meta = read.table('data/spot_meta/clusters_human_mix.txt', sep='\t',
                          header=T, row.names = 1)
# Re-labelling #
spot_starts <- str_sub(rownames(cluster_meta), 1, 2)
mapping <- list("A1"="A1_treated_", "B1"="B1_treated_", 
                "C1"="C1_untreated_", "D1"="D1_untreated_")
spot_ends <- str_sub(rownames(cluster_meta), 4)
spot_names <- c()
for (i in 1:length(spot_starts)) {
  spot_names <- c(spot_names, paste(mapping[[spot_starts[i]]], spot_ends[i], sep=''))
}

rownames(cluster_meta) <- spot_names

cluster_meta <- cluster_meta[names(human$orig.ident),]

# Adding the meta data #
human <- AddMetaData(human, cluster_meta)

################################################################################
                   # Examining the joint UMAP space #
################################################################################
human <- RunPCA(human, assay = "SCT", verbose = FALSE)
human <- RunUMAP(human, reduction = "pca", dims = 1:50)

DimPlot(human, reduction = "umap", group.by = c('tissue_type', 'sample'),
        pt.size=1)

p1 <- DimPlot(human, reduction = "umap", group.by='tissue_type')
p2 <- DimPlot(human, reduction = "umap", group.by='sample')
print(p1+p2)

################################################################################
                    # Running the integration method #
################################################################################
# https://satijalab.org/seurat/articles/integration_introduction.html

# split the dataset into a list of two seurat objects (stim and CTRL)
obj.list <- SplitObject(human, split.by = "sample")

# normalize and identify variable features for each dataset independently
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = obj.list)

# Attempting to remove the human genes ! #
#human_genes <- getSpeciesGenes(human@assays$Spatial@data, 'hg38')
#human_features <- intersect(human_genes, features)
# 2393 human features

# Finding anchors
anchors <- FindIntegrationAnchors(object.list = obj.list)

# this command creates an 'integrated' data assay
human.combined <- IntegrateData(anchorset = anchors)

# Adding in the predicted labels for the cell types 
labels <- c()
probs <- human.combined@assays$human_vladoiu_preds@data
maxs <- probs['max',]
label_set <- setdiff(rownames(probs), c('max'))
probs <- probs[label_set,]
for (celli in 1:ncol(probs)){
  label <- 'unknown'
  if (maxs[celli]>.2) {
    label <- label_set[which(probs[,celli]==maxs[celli])]
  }
  labels <- c(labels, label)
}

human.combined <- AddMetaData(human.combined, maxs, 'vladoiu.scores')
human.combined <- AddMetaData(human.combined, labels, 'vladoiu.labels')

# Now integrating #
DefaultAssay(human.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
human.combined <- ScaleData(human.combined, verbose = FALSE)
human.combined <- RunPCA(human.combined, npcs = 30, verbose = FALSE)
human.combined <- RunUMAP(human.combined, reduction = "pca", dims = 1:30, n.components=2L)
human.combined <- FindNeighbors(human.combined, reduction = "pca", dims = 1:30)
human.combined <- FindClusters(human.combined, resolution = 0.5)

# Checking to see if the clusters better correspond #
DimPlot(human.combined, reduction = "umap", #split.by='sample', 
        group.by = c('vladoiu.labels','seurat_clusters'), pt.size=1, shuffle=T)

p1 <- DimPlot(human.combined, reduction = "umap", group.by='seurat_clusters')
p2 <- DimPlot(human.combined, reduction = "umap", group.by='sample')
print(p1+p2)

DimPlot(human.combined, reduction = "umap", group.by=c('tissue_type', 'seurat_clusters'))

images <- c('slice1', 'slice1_B1_treated.1', 'slice1_C1_untreated.2', 'slice1_D1_untreated.3')
sizes <- c(2.5, 1.5, 1.5, 1.5)
i <- 1
SpatialDimPlot(human.combined, group.by='tissue_type', image=images[1],
               pt.size.factor = 1.5)
SpatialFeaturePlot(human.combined, features='hg38-OLIG2', image=images[4],
               pt.size.factor = 1.5)

## Measure spot-wise similarities to the Vladoiu developing cerebellum dataset. ##
# For performing differential expression after integration, we switch back to the original data
DefaultAssay(human.combined) <- "SCT"
Idents(human.combined) <- 'seurat_clusters'
nk.markers <- FindAllMarkers(human.combined, #grouping.var='treatment',
                                   verbose = FALSE, assay='SCT', only.pos=T,
                                   use.method='wilcox')
unique(nk.markers[nk.markers[,'cluster']==1,'gene'])

cols <- colnames(human.combined@meta.data)
FeaturePlot(human.combined, "hg38-TUBA1B", keep.scale='all',
            pt.size = 1, order=F, split.by='treatment')

DimPlot(human.combined, reduction = "umap", shuffle=T,#split.by='sample', 
        group.by = c('seurat_clusters', 'tissue_type'), pt.size=1, )

DimPlot(human.combined, reduction = "umap", split.by='treatment', 
        group.by = 'vladoiu.labels', pt.size=1, dims=c(1,2))

#human.combined <- readRDS('data/seurat_rds/integrated.rds')
#saveRDS(human.combined, 'data/seurat_rds/integrated.rds')

################################################################################
              # Performing Cell Cycle Scoring #
################################################################################










