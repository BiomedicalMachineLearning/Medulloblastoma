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

out_dir <- 'data/DE_out/Pseudo_TMM_Limma_Voom/'

################################################################################
                    # Loading & formatting data #
################################################################################
seurat_obj <- readRDS('data/seurat_rds/merged.rds')
counts <- seurat_obj@assays$Spatial@counts
#counts <- seurat_obj@assays$SCT@counts
suffix <- 'human_mix_test'
samples <- seurat_obj$sample
species <- seurat_obj$species

# Pseudo-bulking #
sample_names <- unique(samples)
species_names <- c('human', 'mix') #unique(species)
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
cutoff <- .5
drop <- which(apply(cpm(x), 1, max) < cutoff)
x <- x[-drop,] 

detect_bools <- x$counts > 0
detected_genes <- apply(detect_bools, 1, all) # Genes present in each sample 
x <- x[detected_genes,]
dim(x) # number of genes left

# Visualising before and after filtering #
prefix <- paste(out_dir, 'gene_filtering_density', sep='')
plot_name <- paste(prefix, species, suffix, '.pdf', sep='_')
pdf(plot_name)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Unfiltered data", xlab="Log-cpm")
abline(v=cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", sample_names, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
plot <- legend("topright", sample_names, text.col=col, bty="n")
dev.off()

# Showing effect of normalisation #
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

x2 <- x
x2$samples$norm.factors <- 1

#### Density plots form normalisation ####
prefix <- paste(out_dir, 'normalisation_density', sep='')
plot_name <- paste(prefix, species, suffix, '.pdf', sep='_')
pdf(plot_name)
library(RColorBrewer)
nsamples <- ncol(x)
lcpm <- cpm(x2, log=TRUE)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Pre-TMM data", xlab="Log-cpm")
abline(v=cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", sample_names, text.col=col, bty="n")

x2 <- calcNormFactors(x2, method='TMM')  
x2$samples$norm.factors
lcpm <- cpm(x2, log=TRUE)

plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. TMM data", xlab="Log-cpm")
abline(v=cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
plot <- legend("topright", sample_names, text.col=col, bty="n")
dev.off()

# Normalising #
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

x2 <- x
x2$samples$norm.factors <- 1

######## Box plots of RLE values ######
prefix <- paste(out_dir, 'norm_v_unNorm_boxPlot_RLE', sep='')
plot_name <- paste(prefix, species, suffix, '.pdf', sep='_')
pdf(plot_name)
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
rles <- calc_cell_RLE(lcpm)
boxplot(rles, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="RLE")
x2 <- calcNormFactors(x2, method='TMM')  
x2$samples$norm.factors

lcpm <- cpm(x2, log=TRUE)
rles <- calc_cell_RLE(lcpm)
boxplot(rles, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="RLE")
dev.off()

# Normalising #
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

x2 <- x
x2$samples$norm.factors <- 1

#### Density plots for RLE values ####
prefix <- paste(out_dir, 'RLE_densities', sep='')
plot_name <- paste(prefix, species, suffix, '.pdf', sep='_')
pdf(plot_name)
library(RColorBrewer)
nsamples <- ncol(x)
lcpm <- cpm(x2, log=TRUE)
rles <- calc_cell_RLE(lcpm)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
maxy=5
plot(density(rles[,1]), col=col[1], lwd=2, ylim=c(0,maxy), las=2, main="", xlab="")
title(main="A. Pre-TMM data", xlab="RLE")
abline(v=cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(rles[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", sample_names, text.col=col, bty="n")

x2 <- calcNormFactors(x2, method='TMM')  
x2$samples$norm.factors
lcpm <- cpm(x2, log=TRUE)
rles <- calc_cell_RLE(lcpm)

plot(density(rles[,1]), col=col[1], lwd=2, ylim=c(0,maxy), las=2, main="", xlab="")
title(main="B. TMM data", xlab="RLE")
abline(v=cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(rles[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
plot <- legend("topright", sample_names, text.col=col, bty="n")
dev.off()

# Normalising #
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

x2 <- x
x2$samples$norm.factors <- 1

#### Density scatter plots of the RLE values vs mean expression after TMM ####
lcpm <- cpm(x, log=TRUE)
rles <- calc_cell_RLE(lcpm)
orig_names <- colnames(lcpm)
colnames(lcpm) <- str_c(rep('x', ncol(lcpm)), as.character(1:ncol(lcpm)))#str_c(rep('log-cpm-', ncol(lcpm)), orig_names)
colnames(rles) <- str_c(rep('y', ncol(lcpm)), as.character(1:ncol(lcpm)))#str_c(rep('RLE-', ncol(lcpm)), orig_names)
stat_df <- as.data.frame( cbind(lcpm, rles) )

plots <- list()
for (i in 1:length(orig_names)) {
  xi <- paste('x', as.character(i), sep='')
  yi <- paste('y', as.character(i), sep='')
  # form <- paste(y, x, sep='~')
  # h <- hexbinplot(as.formula(form), data=stat_df, colramp=rf,
  #            xlab="log-cpm", ylab="RLE", main=orig_names[i])
  # print(h)
  h <- ggplot(stat_df, aes_string(xi, yi)) + stat_bin2d(bins=50)
  plots[[i]] <- h + ggtitle(orig_names[i]) + xlab('log-cpm') + ylab('RLE')
}

prefix <- paste(out_dir, 'RLEvLogCPM', sep='')
plot_name <- paste(prefix, species, suffix, '.pdf', sep='_')
pdf(plot_name)
plot_grid(plotlist=plots)
dev.off()

## Normalising ##
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

x2 <- x
x2$samples$norm.factors <- 1

######## Box plots of log-cpm values ######
prefix <- paste(out_dir, 'norm_v_unNorm_boxPlot_', sep='')
plot_name <- paste(prefix, species, '.pdf', sep='')
pdf(plot_name)
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2, method='TMM')  
x2$samples$norm.factors

lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")
dev.off()

## Normalising ##
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

x2 <- x
x2$samples$norm.factors <- 1

###### MDS plot ######
prefix <- paste(out_dir, 'norm_v_unNorm_MDS', sep='')
plot_name <- paste(prefix, species, suffix, '.pdf', sep='_')
pdf(plot_name)
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- x$samples$treat
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

plotMDS(lcpm, labels=treat, col=col.group)
title(main="A. Pre-TMM")

lcpm <- cpm(x2, log=TRUE)
col.group <- x$samples$treat
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

plotMDS(lcpm, labels=treat, col=col.group)
title(main="A. TMM")
dev.off()

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

prefix <- paste(out_dir, 'scatter_LFVvsAvgExpr_', sep='')
plot_name <- paste(prefix, species, suffix, '.pdf', sep='_')
pdf(plot_name)
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-1,13))
dev.off()

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










