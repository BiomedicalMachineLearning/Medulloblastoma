# The purpose of this script is to get marker genes for each of the clusters 
# determined from the stlearn SME clustering to help justify the cluster labels
# between the groups are equivalent. 

################################################################################
                    # Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/"
setwd(work_dir)
source('scripts/helpers.R')

library(Seurat)
library(ggplot2)
library(xlsx)

################################################################################
                            # Loading the data #
################################################################################
path_info <- getSpatialPaths()
data_paths <- path_info[[1]]
data_names <- path_info[[2]]

spatials <- readRDSFiles(data_paths, data_names)

# Loading in each of the labels #
meta <- read.table(paste(work_dir, 'data/spot_meta/clusters_human_mix.txt', sep=''),
                   sep='\t', row.names = 1, header=T)

################################################################################
                          # Calling DE genes #
################################################################################
de_by_sample <- list()
for (i in 1:length(spatials)) {
  spatial <- spatials[[i]]
  
  # Subsetting to just data #
  Idents(spatial) <- 'species'
  human <- subset(spatial, idents=c('human', 'mix'))
  human <- SCTransform(human, assay = "Spatial", verbose = T,
                                  return.only.var.genes=F #Ensures all genes returned !
          ) 
  
  
  human_genes <- getSpeciesGenes(human@assays$SCT@data, 'hg38-')
  
  prefix <- paste(str_split(names(spatials)[i], '_')[[1]][1], '-', sep='')
  cell_names <- rownames(human@meta.data)
  prefixed_cell_names <- str_c(rep(prefix, length(cell_names)), cell_names)
  meta_sub <- meta[prefixed_cell_names,]
  rownames(meta_sub) <- cell_names
  human <- AddMetaData(human, meta_sub)
  
  # Calling marker genes with Wilcoxon rank-sum tests #
  de_by_cluster <- list() 
  Idents(human) <- 'tissue_type'
  tissue_types <- unique(human$tissue_type)
  for (j in 1:length(tissue_types)) {
    tissue_type <- tissue_types[j]
    de_df <- FindMarkers(human, test.use='wilcox', features = human_genes,
                                                ident.1 = tissue_type, only.pos=T,
                                                min.pct = .01, logfc.threshold = .05)
    de_bool <- de_df[,'p_val_adj'] < .05
    de_by_cluster[[tissue_type]] <- de_df[de_bool,]
  }
  de_by_sample[[i]] <- de_by_cluster
  
  # Saving the marker genes for each tissue type #
  file_name <- paste('data/DE_out/stlearn_cluster_DE/de_by_cluster.xlsx')
  for (j in 1:length(tissue_types)) {
    sample <- str_split(names(spatials)[i], '_')[[1]][1]
    sheet_name <- paste(sample, tissue_types[j], sep='_')
    sheet_name <- str_replace_all(sheet_name, '/', '-')
    if (sheet_name == "B1_meninges-vasculogenic-mimicry") {
      sheet_name <- 'B1_minges-vasulogenic'
    }
    write.xlsx(de_by_cluster[[j]], file_name,  append=T,
               sheetName=sheet_name, col.names=TRUE, row.names=TRUE)
  }
}

################################################################################
                  # Comparing the DE genes between samples #
################################################################################
library(openxlsx)
# Reading back in the results for comparison #
de_by_sample <- list()
file_name <- paste('data/DE_out/stlearn_cluster_DE/de_by_cluster.xlsx')
sheet_index <- 1
for (i in 1:length(spatials)) {
  prefix <- paste(str_split(names(spatials)[i], '_')[[1]][1], '-', sep='')
  human_bool <- spatials[[i]]$species=='human' |  spatials[[i]]$species=='mix'
  cell_names <- rownames(spatials[[i]]@meta.data)[human_bool]
  prefixed_cell_names <- str_c(rep(prefix, length(cell_names)), cell_names)
  meta_sub <- meta[prefixed_cell_names,]
  
  tissue_types <- unique(meta_sub[,'tissue_type'])
  de_by_cluster <- list()
  for (j in 1:length(tissue_types)) {
    sample <- str_split(names(spatials)[i], '_')[[1]][1]
    sheet_name <- paste(sample, tissue_types[j], sep='_')
    sheet_name <- str_replace_all(sheet_name, '/', '-')
    
    de_df <- read.xlsx(file_name, sheet = sheet_index)
    rownames(de_df) <- de_df[,1]
    de_df <- de_df[,-1]
    
    de_by_cluster[[tissue_types[j]]] <- de_df
    sheet_index <- sheet_index + 1
  }
  
  de_by_sample[[names(spatials)[i]]] <- de_by_cluster
}

sample_gene_overlap <- rownames(spatials[[1]]@assays$SCT@data)
for (i in 2:length(spatials)) {
  sample_gene_overlap <- intersect(sample_gene_overlap, rownames(spatials[[i]]@assays$SCT@data))
}

# Determining overlap of DE genes between the samples #
tissue_types <- unique(meta$tissue_type)
tissue_sims <- list()
tissue_conts <- list()
intersect_scores <- list()
for (i in 1:length(tissue_types)) {
  tissue_type <- str_replace_all(tissue_types[i], '/', '-')
  intersect_scores[[tissue_type]] <- integer(length(sample_gene_overlap))
  names(intersect_scores[[tissue_type]]) <- sample_gene_overlap
  
  samples <- c()
  for (j in 1:length(de_by_sample)) {
    if (tissue_type %in% names(de_by_sample[[j]])) {
       samples <- c(samples, names(de_by_sample)[j])
    }
  }
  
  if (length(samples) == 0) {next}
  
  tissue_conts[[tissue_type]] <- list() 
  sims <- matrix(0, length(samples), length(samples))
  for (l in 1:length(samples)) {
    genes_l <- rownames( de_by_sample[[samples[l]]][[tissue_type]] )
    for (k in 1:length(samples)) {
      genes_k <- rownames( de_by_sample[[samples[k]]][[tissue_type]] )
      
      over_genes <- intersect(genes_l, genes_k)
      diff_genes <- setdiff(genes_l, genes_k)
      overlap <- length( over_genes )
      diff_ <- length( diff_genes )
      
      all_overlap <- 0
      all_diff <- 0
      tissue_intersect <- intersect(names(de_by_sample[[samples[l]]]),
                                    names(de_by_sample[[samples[k]]]))
      for (i_ in 1:length(tissue_intersect)){
        tissue_type_i_ <- tissue_intersect[i_]
        if (tissue_type_i_ != tissue_type) {
          genes_i <- rownames( de_by_sample[[samples[k]]][[tissue_type_i_]] )
          all_overlap <- all_overlap + length( intersect(genes_l, genes_i) )
          all_diff <- all_diff + length( setdiff(genes_l, genes_i) )
        }
      }
      
      cont_table <- matrix(c(overlap, all_overlap, diff_, all_diff),nrow=2,ncol=2)
      colnames(cont_table) <- c('overlap genes', 'diff genes')
      rownames(cont_table) <- c('same cluster', 'other clusters')
      fet_result <- fisher.test(cont_table, alternative="greater")
      sims[l, k] <- fet_result$p.value
      
      cont_name <- paste(samples[l], samples[k], sep="_")
      tissue_conts[[tissue_type]][[cont_name]] <- cont_table
      
      intersect_scores[[tissue_type]][over_genes] <- intersect_scores[[tissue_type]][over_genes] + 1
    }
  }
  
  row.names(sims) <- samples
  colnames(sims) <- samples
  tissue_sims[[tissue_types[i]]] <- sims
}

for (i in 1:length(intersect_scores)) {
  name <- names(intersect_scores[i])
  indices <- order(intersect_scores[[i]], decreasing=T)
  intersect_scores[[name]] <- intersect_scores[[i]][indices]
}





SpatialFeaturePlot(spatial, features = 'hg38-NHLH1')

# TODO perform the DE for each tissue type #





