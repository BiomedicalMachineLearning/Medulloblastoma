### The purpose of this script is to run DE using a Pseudo_TMM_Limma_Voom approach ####
# References:
# * https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# * https://bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#reading-in-count-data
# * https://www.stephaniehicks.com/2018-bioinfosummer-scrnaseq/cleaning-the-expression-matrix.html

################################################################################
                    # Environment setup #
################################################################################
# TODO convert the below to use snakemake inputs #
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/"
setwd(work_dir)
source('scripts/helpers.R')

library(Seurat)
library(ggplot2)
library(edgeR) # Contains limma
library(stringr)
library(hexbin)
library(cowplot)
Sys.setenv(RETICULATE_PYTHON = "/Users/uqbbalde/miniconda3/envs/MedullaBlastoma/bin/python/")
library(reticulate)

out_dir <- 'figure_components/PseudoLimma/DE_figures/'

################################################################################
                    # Loading & formatting data #
################################################################################
seurat_obj <- readRDS('data/seurat_rds/merged.rds')
counts <- seurat_obj@assays$Spatial@counts
#counts <- seurat_obj@assays$SCT@counts
suffix <- 'human'
samples <- seurat_obj$sample
species <- seurat_obj$species

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
            # Normalising data with TMM & Plotting QC plots #
################################################################################
# CPM normalise #
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

# Filtering genes #
#cutoff <- 3.5
#vals <- apply(cpm(x), 1, mean)
#hist(vals[vals<10])
#drop <- which(vals < cutoff)
#x <- x[-drop,]

# Filtering genes #
cutoff <- .5
drop <- which(apply(cpm(x), 1, max) < cutoff)
x <- x[-drop,] 

detect_bools <- x$counts > 0
detected_genes <- apply(detect_bools, 1, all) # Genes present in each sample
x <- x[detected_genes,]
dim(x) # number of genes left

# TMM normalisation #
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

lcpm <- cpm(x, log=TRUE)

#### Density plots form normalisation ####
# RETICULATE NOT WORKING FOR SCIPY #
#source_python('scripts/May_REDO/DE_figure_helpers.py')
#prefix <- paste(out_dir, 'TMM_normalisation_density', sep='')
#plot_name <- paste(prefix, species[1], suffix, '.pdf', sep='_')

#lcpm_vals <- dict( as.list(as.data.frame(lcpm)) )
#colors <- dict( list('A1_treated_human'='firebrick', 'B1_treated_human'='royalblue',
#                     'C1_untreated_human'='limegreen', 'D1_untreated_human'='orchid') )

#density_plot(lcpm_vals, colors, plot_name)

# Saving data for visualisation in python #
lcpm_df <- as.data.frame( lcpm )
write.table(lcpm_df, file='data/DE_out/Pseudo_TMM_Limma_Voom/data/lcpms.txt', sep='\t',
            quote=F)

#### Density plots for RLE values ####
rles <- calc_cell_RLE(lcpm)
rles_df <- as.data.frame( rles )
write.table(rles_df, file='data/DE_out/Pseudo_TMM_Limma_Voom/data/rles.txt', sep='\t',
            quote=F)

################################################################################
            # VOOM mean-variance correction #
################################################################################
design <- model.matrix(~treat)
rownames(design) <- colnames(x$counts)
colnames(design) <- c("Intercept", "treattreated")

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
tfit <- treat(vfit, lfc=0.1)
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

# Saving data necessary to do the scatter plots #
gene_stats <- as.data.frame(tfit$coefficients)
colnames(gene_stats) <- c('logFC')
gene_stats[,'de_status'] <- rep("not-de", nrow(gene_stats))
gene_stats[sig_up,'de_status'] <- 'up'
gene_stats[sig_down,'de_status'] <- 'down'


write.table(gene_stats, file='data/DE_out/Pseudo_TMM_Limma_Voom/data/gene_stats.txt', 
            sep='\t', quote=F)

################ Plotting scatter of LFCvsLogCPMs ##############################
prefix <- paste(out_dir, 'scatter_LFVvsAvgExpr_', sep='')
plot_name <- paste(prefix, species, suffix, '.pdf', sep='_')
pdf(plot_name)
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1],
       xlim=c(-1,13))
dev.off()

#### Density plot of the logFC ############
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
file_name <- paste(out_dir, 'de_results_tfit_', species[1], suffix, '.rds', sep='')
saveRDS(tfit, file=file_name)
#file_name <- paste(out_dir, 'de_results_tfit_human.rds', sep='')
#tfit <- readRDS(file_name)

# Save results #
rank_ <- order(tfit$p.value.adj)
rank_ <- rank_[tfit$p.value.adj[rank_]<cut]
res_df <- as.data.frame(list('t.val'=tfit$t[rank_,1],
                             'log2FoldChange'=tfit$coefficients[rank_,1],
                             'pvalue'=tfit$p.value[rank_,1],
                             'padj'=tfit$p.value.adj[rank_]))

up_df <- res_df[res_df[,'t.val']>0,]
down_df <- res_df[res_df[,'t.val']<0,]

library(xlsx)

file_name <- paste('data/DE_out/Pseudo_TMM_Limma_Voom/de_results_PseudoLimma.xlsx')
sheet_up <- paste(species[1], '.treatment_up', sep='')
sheet_down <- paste(species[1], '.treatment_down', sep='')

write.xlsx(up_df, file_name,
           sheetName=sheet_up, col.names=TRUE, row.names=TRUE, append=T)

write.xlsx(down_df, file_name,
           sheetName=sheet_down, col.names=TRUE, row.names=TRUE, append=T)



