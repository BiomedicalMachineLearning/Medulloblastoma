# The purpose of this script is to perform differential
# expression using DESeq2 between the treated versus untreated
# samples between the different species spots (human, mouse, mix)

# Note that only testing for SCTransform genes, and for each set of spots,
# only testing for the relevant genes.
# human - human genes, mouse - mouse genes, mix - both genes

################################################################################
                    # Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/"
setwd(work_dir)
source('scripts/helpers.R')

library(DESeq2)
library(scran)
library(Seurat)
library(ggplot2)

################################################################################
                        # Loading the data #
################################################################################
# path_info <- getSpatialPaths()
# data_paths <- path_info[[1]]
# data_names <- path_info[[2]]
# 
# spatials <- readRDSFiles(data_paths, data_names)

merged <- readRDS('data/seurat_rds/merged.rds')

################################################################################
              # Calling DE genes on the merged data #
################################################################################
# NOTE: get very stuck on this, found that I needed to use a nested design as
#     shown here:
#     https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#group-specific-condition-effects-individuals-nested-within-groups
# And to get this to work, I needed to specify replicate 1/2 to get it to work.

# -------------------- Method 1: nested design matrix -----------------------#

########## Formatting the data ##########
sub <- merged
print(unique(sub$species))
sctransform_genes <- rownames(sub@assays$SCT@data)
count_matrix <- sub@assays$Spatial@counts[sctransform_genes,]
design_matrix <- as.data.frame(cbind(sub$species, sub$sample, sub$treatment))
colnames(design_matrix) <- c('species', 'sample', 'treatment')
print(design_matrix)

# Labelling replicates #
samples <- design_matrix[,'sample']
replicates <- character(ncol(count_matrix))
replicates[samples=='A1_treated'] <- '1'
replicates[samples=='B1_treated'] <- '2'
replicates[samples=='C1_untreated'] <- '1'
replicates[samples=='D1_untreated'] <- '2'

design_matrix <- DataFrame(species=relevel(factor(design_matrix[,'species']), 'mix'),
                           sample=factor(design_matrix[,'sample']),
                           treatment=relevel(factor(design_matrix[,'treatment']), 'untreated'),
                           replicate=factor(replicates)
)
m1 <- model.matrix(~species + species:replicate + species:treatment, design_matrix)
all.zero <- apply(m1, 2, function(x) all(x==0))
print(all.zero)

########## Calling differential expression ###########
# Following recommendations from DESeq2 vignette undersection:
# 'Recommendations for single cell analysis'
count_mat <- count_matrix+1
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = design_matrix,
                              design = ~species + species:replicate + species:treatment)

dds <- DESeq(dds, useT=T, minmu=1e-6, minReplicatesForReplace=Inf)

########### Getting the outputted DEResults and formatting ############
source('scripts/species_DE_treatment_helpers.R')

avail_contrasts <- resultsNames(dds)
print(avail_contrasts)

species_treatments <- c('speciesmouse.treatmenttreated', 
                        'specieshuman.treatmenttreated',
                        'speciesmix.treatmenttreated')
all_de_results <- list()
specFiltered_de_results <- list() # DE where filter genes based on species #
for (i in 1:length(species_treatments)) {
  spec_treat <- species_treatments[i]
  de_results <- results(dds, contrast=list(spec_treat))
  
  # Formatting all the results, regardless of species #
  all_de_result <- formatResults(de_results, padj_cutoff=.05)
  all_de_results[[ spec_treat ]] <- all_de_result
}

# Writing out results to excel #
library(xlsx)

file = 'data/supps/de_results_MDB.xlsx'
for (i in 1:length(all_de_results)) {
  de <- all_de_results[[i]]
  spec_treat <- species_treatments[i]
  spec_treat <- trimws(spec_treat, whitespace = 'species')
  spec_treat <- trimws(spec_treat, whitespace = 'treated')
  for (j in 1:2) {
    sheet_name <- paste(spec_treat, '_', names(de)[j], sep='')
    print(sheet_name)
    if (nrow(de[[j]]) != 0) {
      write.xlsx(de[[j]], file, sheetName=sheet_name, 
               col.names=TRUE, row.names=TRUE, append=T)
    }
  }
}

saveRDS(dds, 'data/DESeq_object.rds')



