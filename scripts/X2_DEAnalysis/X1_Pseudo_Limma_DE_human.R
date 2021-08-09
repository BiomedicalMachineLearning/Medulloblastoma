### The purpose of this script is to run DE using a Pseudo_TMM_Limma_Voom approach ####
# References:
# * https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# * https://bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#reading-in-count-data
# * https://www.stephaniehicks.com/2018-bioinfosummer-scrnaseq/cleaning-the-expression-matrix.html

################################################################################
                    # Environment setup #
################################################################################

library(Seurat)
library(ggplot2)
library(edgeR) # Contains limma
library(stringr)
library(hexbin)
library(cowplot)
library(reticulate)
pd <- import('pandas')

work_dir <- rstudioapi::getSourceEditorContext()$path
work_dir <- str_split(work_dir, 'scripts/')[[1]][1]
setwd(work_dir)

source('scripts/utils/helpers.R')
source('scripts/X2_DEAnalysis/limma_helpers.R')

suffix <- 'human'
out_dir <- 'data/DE_out/Pseudo_Limma_human/'
species_dir <- 'data/spot_meta/species_classify_v2/'
stats_dir <- paste(out_dir, 'data/', sep='')

################################################################################
                    # Loading & formatting data #
################################################################################
seurat_obj <- readRDS('data/seurat_rds/all.rds')
counts <- seurat_obj@assays$Spatial@counts

# Loading species from file #
samples <- c('A1', 'B1', 'C1', 'D1')
prefixes <- c('A1_treated_', 'B1_treated_', 'C1_untreated_', 'D1_untreated_')
species2 <- character()
for (i in 1:length(samples)) {
  df <- pd$read_csv(paste0(species_dir, samples[i], '_species.txt'), 
                    sep='\t', index_col=as.integer(0))
  rownames(df) <- str_c(rep(prefixes[i], nrow(df)), rownames(df))
  old_names <- names(species2)
  species2 <- c(species2, df[,'species'])
  names(species2) <- c(old_names, rownames(df))
}
print(all(names(seurat_obj$orig.ident)==names(species2)))

species <- species2
samples <- seurat_obj$sample

# Pseudo-bulking #
sample_names <- unique(samples)
species_names <- c('human') #unique(species)
pseudobulk_counts <- list()
for (i in 1:length(sample_names)) {
  sample_name <- sample_names[i]
  sample_bool <- samples==sample_name
  for (j in 1:length(species_names)) {
    species_bool <- species == species_names[j] 
    spot_bool <- sample_bool & species_bool
    pseudobulk <- apply(counts[,spot_bool], 1, sum)
    expr_name <- paste(sample_name, species_names[j], sep='_')
    pseudobulk_counts[[expr_name]] <- pseudobulk
  }
}
pseudobulk_counts <- as.data.frame(pseudobulk_counts)

x <- DGEList(pseudobulk_counts)

# Adding in sample meta data #
split_df <- as.data.frame(str_split(colnames(pseudobulk_counts), '_'))
treat <- as.factor( as.character(split_df[2,]) )
treat <- relevel(treat, 'untreated')
species <- as.factor( as.character(split_df[3,]) )
x$samples$treat <- treat
x$samples$species <- species

################################################################################
  # Saving lcpms prior to filtering so can perform PCA on all samples #
################################################################################
x2 <- x
x2$samples$norm.factors <- 1

x2 <- calcNormFactors(x2, method='TMM')  
x2$samples$norm.factors
lcpm <- cpm(x2, log=TRUE)

# Saving data for visualisation in python #
lcpm_df <- as.data.frame( lcpm )
write.table(lcpm_df, file=paste(stats_dir, 'lcpms_all.txt', sep=''),
            sep='\t', quote=F)

################################################################################
            # Normalising data with TMM & Plotting QC plots #
################################################################################
x <- DGEList(pseudobulk_counts)

# Adding in sample meta data #
split_df <- as.data.frame(str_split(colnames(pseudobulk_counts), '_'))
treat <- as.factor( as.character(split_df[2,]) )
treat <- relevel(treat, 'untreated')
species <- as.factor( as.character(split_df[3,]) )
x$samples$treat <- treat
x$samples$species <- species

# CPM normalise #
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE, prior.count=2)

# Filtering genes #
cutoff <- 2
drop <- which(apply(cpm(x), 1, max) < cutoff)
x <- x[-drop,] 

detect_bools <- x$counts > 0
detected_genes <- apply(detect_bools, 1, all) # Genes present in each sample 
human_genes <- str_detect(rownames(x), 'hg38-')
names(human_genes) <- rownames(x)
genes_bool <- detected_genes & human_genes
x <- x[genes_bool,]
dim(x) # number of genes left
# 13069 genes remaining.

##### Visualising before and after filtering #####
plot_filtering_density(lcpm, x, out_dir, species, suffix)

##### Density plots from normalisation #####
plot_norm_density(x, out_dir, species, suffix)

######## Box plots of RLE values ######
plot_RLE_boxes(x, out_dir, species, suffix)

#### Density plots for RLE values ####
plot_RLE_densities(x, out_dir, species, suffix)

################################################################################
          # Outputting data for plotting results in python #
################################################################################
# TMM normalisation #
x <- calcNormFactors(x, method = "TMM")
lcpm <- cpm(x, log=TRUE)

# Saving data for visualisation in python #
lcpm_df <- as.data.frame( lcpm )
write.table(lcpm_df, file=paste(stats_dir, 'lcpms.txt', sep=''),
            sep='\t', quote=F)

#### Density plots for RLE values ####
rles <- calc_cell_RLE(lcpm)
rles_df <- as.data.frame( rles )
write.table(rles_df, file=paste(stats_dir, 'rles.txt', sep='') ,
            sep='\t', quote=F)

################################################################################
            # VOOM mean-variance correction #
################################################################################
design <- model.matrix(~treat)
rownames(design) <- colnames(x$counts)

contr.matrix <- makeContrasts(
  treatment = treattreated, 
  levels = colnames(design))
contr.matrix

#### Plotting the correction of the mean-variance tread ####
prefix <- paste(out_dir, 'voom_meanvar_correction', sep='')
plot_name <- paste(prefix, species, suffix, '.pdf', sep='_')
pdf(plot_name)
par(mfrow=c(1,2))
x <- calcNormFactors(x, method = "TMM")
v <- voom(x, design, plot=TRUE)

# Performing the mean-variance correction #
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)

plotSA(efit, main="Final model: Mean-variance trend")
dev.off()

################################################################################
                # Calling Differential Expression #
################################################################################
# Without the TREAT criterion #
summary(decideTests(efit))

# With the TREAT criterion #
tfit <- treat(vfit, lfc=0.15)
dt <- decideTests(tfit)

# Getting significant genes #
tfit$p.value.adj <- p.adjust(tfit$p.value, method='fdr')
dt_adj <- dt

cut <- .05
sig_up <- tfit$p.value.adj < cut & tfit$t[,1] > 0
sig_down <- tfit$p.value.adj < cut & tfit$t[,1] < 0
not_sig <- tfit$p.value.adj > cut

dt_adj[sig_up,1] <- 1
dt_adj[sig_down,1] <- -1
dt_adj[not_sig,1] <- 0
summary(dt_adj)
# treatment
# Down         511
# NotSig     11869
# Up           689

prefix <- paste(out_dir, 'scatter_LFVvsAvgExpr_', sep='')
plot_name <- paste(prefix, species, suffix, '.pdf', sep='_')
pdf(plot_name)
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-1,13))
dev.off()

####### Saving data necessary to do the scatter plots #######
gene_stats <- as.data.frame(tfit$coefficients)
colnames(gene_stats) <- c('logFC')
gene_stats[,'de_status'] <- rep("not-de", nrow(gene_stats))
gene_stats[sig_up,'de_status'] <- 'up'
gene_stats[sig_down,'de_status'] <- 'down'

write.table(gene_stats, file=paste(stats_dir, 'gene_stats.txt', sep=''),
            sep='\t', quote=F)

#### Density plot of the logFC 
prefix <- paste(out_dir, '_density', sep='')
plot_name <- paste(prefix, 'logFC', suffix, '.pdf', sep='_')
pdf(plot_name)
library(RColorBrewer)
nsamples <- ncol(x)
lcpm <- cpm(x2, log=TRUE)
col <- brewer.pal(3, "Paired")
par(mfrow=c(1,1))
plot(density(tfit$coefficients[,1]), col=col[1], lwd=2, ylim=c(0,1), las=2, main="", xlab="")
title(main="logFC density", xlab="logFC")
dev.off()

##### Saving the results #####
file_name <- paste(out_dir, 'de_results_tfit_', '', suffix, '.rds', sep='')
saveRDS(tfit, file=file_name)
#file_name <- paste(out_dir, 'de_results_tfit_human.rds', sep='')
#tfit <- readRDS(file_name)

# Save results #
# Function is from limma_helpers.R #
# tfit <- readRDS(file_name)
save_results(tfit, out_dir, 'human', suffix, save=T)









