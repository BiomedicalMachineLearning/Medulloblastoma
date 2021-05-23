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
source('scripts/conditionSpatialDE/scripts/limma_helpers.R')

library(Seurat)
library(ggplot2)
library(edgeR) # Contains limma
library(stringr)
library(hexbin)
library(cowplot)

out_dir <- 'data/DE_out/Spot_TMM_Limma_Voom/'

################################################################################
                    # Loading & formatting data #
################################################################################
seurat_obj <- readRDS('data/seurat_rds/merged.rds')
counts <- seurat_obj@assays$SCT@counts
#counts <- seurat_obj@assays$SCT@counts
suffix <- ''
samples <- seurat_obj$sample
species_spot <- seurat_obj$species

species <- 'human'
spec_bool <- species_spot==species
counts <- counts[,spec_bool]
x <- DGEList(counts)

# Adding in sample meta data #
treat <- as.factor( seurat_obj$treatment[spec_bool] )
treat <- relevel(treat, 'untreated')
x$samples$treat <- treat
samples <- samples[spec_bool]
sample_names <- unique(samples)

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
# Genes present in each sample 
detected_genes <- apply(detect_bools, 1, 
                        function(x) {for (i in 1:length(sample_names)) {
                                        if (!any(x[samples==sample_names[i]])) {return (F)}
                                      }; return (T)})
x <- x[detected_genes,]
dim(x) # number of genes left

# Normalising #
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

x2 <- x
x2$samples$norm.factors <- 1

######## Box plots of RLE values ######
prefix <- paste(out_dir, 'norm_v_unNorm_boxPlot_RLE', sep='')
plot_name <- paste(prefix, species, suffix, '.pdf', sep='_')
pdf(plot_name)

sample_rles <- calcSampleRLE(x2, samples, sample_names, mean)

par(mfrow=c(1,2))
boxplot(sample_rles, las=2, col=col, main="")
title(main="A. Unnormalised data",ylab="RLE")

x2 <- calcNormFactors(x2, method='TMM')  

sample_rles <- calcSampleRLE(x2, samples, sample_names, mean)

boxplot(sample_rles, las=2, col=col, main="")
title(main="B. TMM data",ylab="RLE")

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

dt_adj[sig_up,1] <- 1
dt_adj[sig_down,1] <- -1
summary(dt_adj)

prefix <- paste(out_dir, 'scatter_LFVvsAvgExpr', sep='')
plot_name <- paste(prefix, species, suffix, '.pdf', sep='_')
pdf(plot_name)
plotMD(tfit, column=1, status=dt_adj[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))
dev.off()

##### Saving the results #####
saveRDS(tfit, file="data/DE_out/TMM_Limma_Voom/de_results_tfit.rds")

# Save results #
rank_ <- order(tfit$p.value.adj)
rank_ <- rank_[tfit$p.value.adj[rank_]<cut]
res_df <- as.data.frame(list('t.val'=tfit$t[rank_,1], 
                             'pvalue'=tfit$p.value[rank_,1],
                             'padj'=tfit$p.value.adj[rank_]))

up_df <- res_df[res_df[,'t.val']>0,]
down_df <- res_df[res_df[,'t.val']<0,]

library(xlsx)

file_name <- paste('data/DE_out/Spot_TMM_Limma_Voom/de_results_SpotLimma.xlsx')
sheet_up <- paste(species[1], '.treatment_up', sep='')
sheet_down <- paste(species[1], '.treatment_down', sep='')

write.xlsx(up_df, file_name, 
           sheetName=sheet_up, col.names=TRUE, row.names=TRUE, append=T)

write.xlsx(down_df, file_name, 
           sheetName=sheet_down, col.names=TRUE, row.names=TRUE, append=T)










