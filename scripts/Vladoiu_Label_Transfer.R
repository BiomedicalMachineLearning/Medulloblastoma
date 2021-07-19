# Performs the label transfer using the 
# Vladoiu as a reference, and the spatial spots
# as the target !
#   INPUT: -> data/third_party_data/Vladoiu2019_Nature_scRNA/vladoiu.rds,
#          -> data/seurat_rds/all.rds
# OUTPUT:  -> data/seurat_rds/all.rds

################################################################################
                              # Environment setup #
################################################################################

library(Seurat)
library(ggplot2)
source('scripts/MedullaBlastoma/helpers.R')

data_dir <- '/30days/uqbbalde/MedullaBlastoma/'

################################################################################
                              # Loading the data #
################################################################################
#### The spatial data #####
spatial_dir <- paste(data_dir, 'data/seurat_rds/all.rds', sep='')
merged <- readRDS(spatial_dir)

### Separating into data v mouse genes in order to get different
### spaces, that way can represent the variation due to mouse cells
### versus variation due to data expression. Will give different
### interpretations. 
human_gene_names <- getSpeciesGenes(merged@assays$Spatial@counts, 'hg38-') 
mouse_gene_names <- getSpeciesGenes(merged@assays$Spatial@counts, 'mm10-')

human_counts <- merged@assays$Spatial@counts[human_gene_names,]
rownames(human_counts) <- convertSpeciesGenes(human_gene_names)

mouse_counts <- merged@assays$Spatial@counts[mouse_gene_names,]
rownames(mouse_counts) <- str_replace(mouse_gene_names, 'mm10-', "")

human <- CreateSeuratObject(human_counts)
mouse <- CreateSeuratObject(mouse_counts) 

# Now performing the SCTransform normalisation & dim reduction #
human <- SCTransform(human, verbose=T, return.only.var.genes=F) %>% 
	 RunPCA() %>% RunUMAP(dims = 1:30)
mouse <- SCTransform(mouse, verbose=T, return.only.var.genes=F) %>%
         RunPCA() %>% RunUMAP(dims = 1:30)

#### The scRNA-seq data #####
mouse_dir <- paste(data_dir, 'Vladoiu2019_Nature_scRNA/vladoiu2019.rds', sep='')
mouse_dev <- readRDS(mouse_dir)

################################################################################
                   # Performing the label transfer #
################################################################################
########### Human ###########
## Calculating the anchors ##
human_anchors <- FindTransferAnchors(reference = mouse_dev, 
				query = human, normalization.method = "SCT")

## Making cell type predictions ##
human_preds.assay <- TransferData(anchorset = human_anchors, 
				refdata = mouse_dev$cell_labels, prediction.assay = TRUE, 
    				weight.reduction = human[["pca"]], 
				dims = 1:30)
merged[['human_vladoiu_preds']] <- human_preds.assay

########### Mouse ###########
## Calculating the anchors ##
mouse_anchors <- FindTransferAnchors(reference = mouse_dev,
                                query = mouse, normalization.method = "SCT")

## Making cell type predictions ##
mouse_preds.assay <- TransferData(anchorset = mouse_anchors,
                                refdata = mouse_dev$cell_labels, prediction.assay = TRUE,
                                weight.reduction = mouse[["pca"]],
                                dims = 1:30)
merged[['mouse_vladoiu_preds']] <- mouse_preds.assay

################################################################################
             # Saving the results for download & visualisation #
################################################################################
saveRDS(merged, spatial_dir)

# Saving the text files containing the cell estimates #
decon_dir <- paste(data_dir, 'data/decon_out/', sep='')

# Deconvolution results based on using the data gene homologues #
human_path <- paste(decon_dir, 'human_vladoiu_cellprops.txt', sep='')
write.table(merged@assays$human_vladoiu_preds@data, human_path, 
		sep='\t', quote=F)

# Deconvolution results based on using mouse gene homologues #
mouse_path <- paste(decon_dir, 'mouse_vladoiu_cellprops.txt', sep='')
write.table(merged@assays$mouse_vladoiu_preds@data, mouse_path,
                sep='\t', quote=F)







