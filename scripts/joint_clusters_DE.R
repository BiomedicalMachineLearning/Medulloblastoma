# Purpose of this script is to call DE genes for the clusters via pseudobulking #

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

########################################################################################
# TODO measure spot-wise similarities to the Vladoiu developing cerebellum dataset.
# TODO perform NNMF to create factors representing combinations of these cell types.
# TODO plot the factor to see if the clusters represent different combinations of cell types.
########################################################################################
human.combined <- readRDS('data/seurat_rds/integrated.rds')

spot_nmf <- read.csv('data/spot_meta/vladoiu_spot_nmf_scores.txt',
                     sep='\t', row.names=1)

spot_scores <- read.csv('data/spot_meta/vladoiu_spot_scores.txt', sep='\t',
                        row.names=1)

human.combined <- AddMetaData(human.combined, spot_nmf)

human.combined <- AddMetaData(human.combined, spot_scores)

DefaultAssay(human.combined) <- "SCT"
DimPlot(human.combined, reduction = "umap", group.by='seurat_clusters')

cells <- colnames(spot_scores)
FeaturePlot(human.combined, 'hg38-AL513365.2', pt.size = 1, order=T)

############################################################################################
# TO get DE genes, will try pseudobulking of cells of all the samples & compare to rest... #
############################################################################################
counts <- human.combined@assays$SCT@counts

samples <- human.combined$sample
clusters <- human.combined$seurat_clusters

# Pseudo-bulking #
sample_names <- unique(samples)
cluster_names <- unique(clusters)
pseudobulk_counts <- list()
for (i in 1:length(sample_names)) {
  sample_name <- sample_names[i]
  sample_bool <- samples==sample_name
  for (j in 1:length(cluster_names)) {
    cluster_bool <- clusters == cluster_names[j] 
    spot_bool <- sample_bool & cluster_bool
    pseudobulk <- apply(counts[,spot_bool], 1, sum)
    expr_name <- paste(sample_name, cluster_names[j], sep='_')
    pseudobulk_counts[[expr_name]] <- pseudobulk
  }
}
pseudobulk_counts <- as.data.frame(pseudobulk_counts)

# Pseudo-bulking for one cluster versus the rest #
sample_names <- unique(samples)
cluster_name <- 0
pseudobulk_counts <- list()
for (i in 1:length(sample_names)) {
  sample_name <- sample_names[i]
  sample_bool <- samples==sample_name
  # Cluster #
  cluster_bool <- clusters == cluster_name
  spot_bool <- sample_bool & cluster_bool
  pseudobulk <- apply(counts[,spot_bool], 1, sum)
  expr_name <- paste(sample_name, cluster_name, sep='_')
  pseudobulk_counts[[expr_name]] <- pseudobulk
  # Rest #
  rest_bool <- clusters != cluster_name
  spot_bool <- sample_bool & rest_bool
  pseudobulk <- apply(counts[,spot_bool], 1, sum)
  expr_name <- paste(sample_name, 'rest', sep='_')
  pseudobulk_counts[[expr_name]] <- pseudobulk
}
pseudobulk_counts <- as.data.frame(pseudobulk_counts)

x <- DGEList(pseudobulk_counts)

# Adding in sample meta data #
split_df <- as.data.frame(str_split(colnames(pseudobulk_counts), '_'))
sample <- as.factor( as.character(split_df[2,]) )
sample <- relevel(sample, 'untreated')
cluster <- as.factor( as.character(split_df[3,]) )
x$samples$sample <- sample
x$samples$cluster <- cluster

######################## Gene filtering ########################
# CPM normalise #
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

# Filtering genes #
cutoff <- 3
drop <- which(apply(cpm(x), 1, max) < cutoff)
x <- x[-drop,] 
dim(x)

min_samples <- function(vals) {sum(vals)>3}

detect_bools <- x$counts > 15
detected_genes <- apply(detect_bools, 1, min_samples) # Genes present in each sample 
x <- x[detected_genes,]
dim(x) # number of genes left

# Checking the logcpm distributions #
library(RColorBrewer)
nsamples <- ncol(x)
lcpm <- cpm(x, log=TRUE)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Pre-TMM data", xlab="Log-cpm")
abline(v=cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(x), text.col=col, bty="n")

################################################################################
              # VOOM mean-variance correction #
################################################################################
design <- model.matrix(~cluster)
rownames(design) <- colnames(x$counts)

contr.matrix <- makeContrasts(
  cluster1 = cluster3, 
  levels = colnames(design))
contr.matrix

#### Plotting the correction of the mean-variance tread ####
v <- voom(x, design, plot=TRUE)

# Performing the mean-variance correction #
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)

plotSA(efit, main="Final model: Mean-variance trend")

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

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-1,13))

# Looking at results #
rank_ <- order(tfit$p.value.adj)
rank_ <- rank_[tfit$p.value.adj[rank_]<cut]
res_df <- as.data.frame(list('t.val'=tfit$t[rank_,1], 
                             'log2FoldChange'=tfit$coefficients[rank_,1],
                             'pvalue'=tfit$p.value[rank_,1],
                             'padj'=tfit$p.value.adj[rank_]))

up_df <- res_df[res_df[,'t.val']>0,]
down_df <- res_df[res_df[,'t.val']<0,]


