####### From calling differential expression ! ##########



# --- Method 1: within species, DE between treat v control, sample agnostic ---#
# Idents(merged) <- 'species'
# species <- 'data'
# sub <- subset(merged, idents = species)

######### From incorrect attempt ####
# species_names <- c('data', 'mouse', 'mix')
# sample_names <- c('A1_treated', 'B1_treated', 'C1_untreated', 'D1_untreated')
# sample_nested <- character(ncol(count_matrix))
# for (i in 1:length(species_names)) {
#   species_bool <- design_matrix[,'species'] == species_names[i]
#   
#   for (j in 1:length(sample_names)) {
#     sample_bool <- design_matrix[,'sample'] == sample_names[j]
#     sample_nested[species_bool&sample_bool] <- as.character(i)
#   }
# }

######## Data formatting ########
# sub <- merged
# sctransform_genes <- rownames(sub@assays$SCT@data)
# count_matrix <- sub@assays$Spatial@counts[sctransform_genes,]
# design_matrix <- as.data.frame(cbind(sub$sample, sub$treatment, sub$species))
# colnames(design_matrix) <- c('sample', 'treatment', 'species')
# 
# # Resolving issue 'Matrix not full rank' by using a nested design,
# # as described here: 
# # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#model-matrix-not-full-rank
# # WHERE:
# # sample = individual, species = group, treatment = condition
# design_matrix <- DataFrame(species=factor(design_matrix[,'species']),
#                            sample=factor(design_matrix[,'sample']),
#                            treatment=factor(design_matrix[,'treatment'])
#                           )
# model_mat <- model.matrix(~treatment + treatment:species + treatment:sample, coldata)
# 
# treatment <- as.character(design_matrix[,'treatment']) == "treated"
# mouse <- as.character(design_matrix[,'species']) == "mouse"
# data <- as.character(design_matrix[,'species']) == "data"
# 
# samples_nested <- character( length(treatment) )
# samples_nested <- 
# 
# design_matrix$sample.nested <- factor() 
# 
# dds <- DESeqDataSetFromMatrix(countData = count_matrix,
#                               colData = design_matrix,
#                               design = ~ species + species:sample + species:treatment)
# 

# From when I was going to limit the DEGenes to just the #
# genes that were in a particular species when calling DE #

# Idents(merged) <- 'species'
# species <- 'data'
# sub <- subset(merged, idents = species) 
# print(unique(sub$species))
# sctransform_genes <- rownames(sub@assays$SCT@data)
# count_matrix <- sub@assays$Spatial@counts[sctransform_genes,]
# design_matrix <- as.data.frame(cbind(sub$sample, sub$treatment))
# colnames(design_matrix) <- c('sample', 'treatment')
# 
# if (species == 'data') {
#   sub_genes <- getSpeciesGenes(count_matrix, 'hg38')
# } else if (species == 'mouse') { 
#   sub_genes <- getSpeciesGenes(count_matrix, 'mm10')
# } else {
#   sub_genes <- sctransform_genes
# }
# 
# count_matrix <- count_matrix[sub_genes,]
# 






####### From predicting data vs mouse ##########

#### NOT IN USE: Testing with on case ####
# human_scores <- getSpeciesScore(spatials[[1]], species_prefix='hg38', 
#                                 assay='SCT', matrix='data') 
# mouse_scores <- getSpeciesScore(spatials[[1]], species_prefix='mm10',
#                                 assay='SCT', matrix='data') 
# 
# scores = as.data.frame(cbind(human_scores, mouse_scores))
# 
# # Setting cutoffs
# human_cutoffs <- c(250, 400)
# mouse_cutoffs <- c(250, 250)
# species <- character(nrow(scores))
# 
# human_bool_x <- scores[,1]>human_cutoffs[1]
# human_bool_y <- scores[,2]<human_cutoffs[2]
# species[human_bool_x & human_bool_y] <- 'data'
# 
# mouse_bool_x <- scores[,1]<mouse_cutoffs[1]
# mouse_bool_y <- scores[,2]>mouse_cutoffs[2]
# species[mouse_bool_x & mouse_bool_y] <- 'mouse'
# 
# species[species==""] <- 'mix'
# 
# scores[,'species'] <- species
# ggplot(scores, aes(x=mouse_scores, y=human_scores, color=species)) + geom_point()

###################### Junk code ###########################
# Plotting the scores #
# Total reads per spot indicate that data requires further filtering.
# Will need to get the spots used by Ryan. 
# In the meantime, just going
# to filter to generate a putative list. 
#ggplot(scores, aes(human_scores+mouse_scores)) + geom_histogram() 


#scores[,'clusters'] <- as.character(clusters)
#mapping <- list('0'='data', '1'='mouse', '2'='mouse',
#                '3'='mouse', '4'='data', '5'='mix')
#scores[,'species'] <- replaceByMapping(mapping, scores[,'clusters'])


#clusterer <- sklearn_mixture$GaussianMixture(n_components=as.integer(2), 
#                                             covariance_type='spherical')
#train_indices <- sample.int(nrow(scores), nrow(scores)/5, replace=F)
#clusterer$fit(scores[train_indices, 1:2])
#probs <- clusterer$predict_proba( scores[,1:2] )[,1] # NOTE only using probability of one cluster
#scores[,'probs'] <- probs
#ggplot(scores, aes(probs)) + geom_histogram() 
# Using .05 cutoff to determine which below to which cluster #
#clusters <- character(length(probs))
#cutoff <- .01
#clusters[probs>1-cutoff] <- '1'
#clusters[probs<cutoff] <- '2'
#clusters[clusters==""] <- 'mix'

# Setting cutoffs #
#ratios <- mouse_scores/human_scores
#ratios[is.infinite(ratios)] <- max(ratios[is.infinite(ratios)==F])
#hist(ratios[ratios>0.2], breaks=100)
#scores[,'ratios'] <- ratios

# Clustering into groups to identify data vs mouse vs mixture #
#clusterer <- sklearn_cluster$KMeans(n_clusters=as.integer(2))
#clusters <- clusterer$fit_predict(scores[,1:2])
#scores[,'clusters'] <- clusters
# ggplot(scores, aes(probs, color=clusters)) + geom_histogram() 
# print(length(unique(clusters)))
# 
# classifier <- sklearn_nbayes$GaussianNB()
# classifier$fit(scores[,1:2], scores[,'clusters'])
# probs <- classifier$predict_proba(scores[,1:2])[,1]
# scores[,'probs'] <- probs
# clusters <- character(length(probs))
# cutoff <- .01
# clusters[probs>1-cutoff] <- '1'
# clusters[probs<cutoff] <- '2'
# clusters[clusters==""] <- 'mix'
# scores[,'clusters'] <- clusters

#### Adding as new labels to the data ######







