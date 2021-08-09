# Purpose of this script is to impliment helpers for performing the Limma DE #

library(xlsx)
library(RColorBrewer)
cutoff <- .5

calc_cell_RLE <-
  function (expr_mat, spikes = NULL)
    # Calculating the RLE values #
  {
    RLE_gene <- function(x) {
      #if (median(unlist(x)) > 0) {
        log((x + 1)/(median(unlist(x)) + 1))/log(2)
      #}
      #else {
      #  rep(NA, times = length(x))
      #}
    }
    
    if (!is.null(spikes)) {
      RLE_matrix <- t(apply(expr_mat[-spikes, ], 1, RLE_gene))
    }
    else {
      RLE_matrix <- t(apply(expr_mat, 1, RLE_gene))
    }
    #cell_RLE <- apply(RLE_matrix, 2, median, na.rm = T)
    return(RLE_matrix)
  }

perSampleRLE <- function(spot_rles, samples, sample_names) {
  # Calculates the per sample RLE by subsamplying so each has the
  # same number of spots and getting the mean RLE for each spot
  getSampleSizes <- function(samples){ 
    sample_sizes <- as.integer(length(sample_names))
    for (i in 1:length(sample_names)) {
      sample_sizes[i] <- sum(samples==sample_names[i])
    }
    return (sample_sizes)
  }
  min_sample_size <- min(getSampleSizes(samples))
  
  sample_rles <- list()
  for (i in 1:length(sample_names)){
    indices <- sample(x=which(samples==sample_names[i]), size=min_sample_size, replace=F)
    sample_rles[[sample_names[i]]] <- spot_rles[indices] 
  }
  sample_rles <- as.data.frame(sample_rles)
}

calcSampleRLE <- function(x, samples, sample_names, function_){
  # function_ refers to either mean or median for getting per spot RLE value #
  lcpm <- cpm(x, log=TRUE)
  rles <- calc_cell_RLE(lcpm)
  spot_rles <- apply(rles, 2, median)
  sample_rles <- perSampleRLE(spot_rles, samples, sample_names)
  return(sample_rles)
}

plot_filtering_density <- function(lcpm, x, out_dir, species, suffix){
  # Plots the density before and after filtering #
  # Visualising before and after filtering #
  cutoff <- .5
  prefix <- paste(out_dir, 'gene_filtering_density', sep='')
  plot_name <- paste(prefix, species, suffix, '.pdf', sep='_')
  pdf(plot_name)
  library(RColorBrewer)
  nsamples <- ncol(x)
  col <- brewer.pal(nsamples, "Paired")
  par(mfrow=c(1,2))
  ymax <- .5
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,ymax), las=2, main="", xlab="")
  title(main="A. Unfiltered data", xlab="Log-cpm")
  abline(v=cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", sample_names, text.col=col, bty="n")
  lcpm <- cpm(x, log=TRUE)
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,ymax), las=2, main="", xlab="")
  title(main="B. Filtered data", xlab="Log-cpm")
  abline(v=cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  plot <- legend("topright", sample_names, text.col=col, bty="n")
  dev.off()
  print("Outputted to:")
  print(plot_name)
}

plot_norm_density <- function(x, out_dir, species, suffix) {
  # Plots the normalisation density & saves to the desired directory #
  # Showing effect of normalisation #
  x <- calcNormFactors(x, method = "TMM")
  x$samples$norm.factors
  
  x2 <- x
  x2$samples$norm.factors <- 1
  
  ymax <- .5
  prefix <- paste(out_dir, 'normalisation_density', sep='')
  plot_name <- paste(prefix, species, suffix, '.pdf', sep='_')
  pdf(plot_name)
  nsamples <- ncol(x)
  lcpm <- cpm(x2, log=TRUE)
  col <- brewer.pal(nsamples, "Paired")
  par(mfrow=c(1,2))
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,ymax), las=2, main="", xlab="")
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
  
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,ymax), las=2, main="", xlab="")
  title(main="B. TMM data", xlab="Log-cpm")
  abline(v=cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  plot <- legend("topright", sample_names, text.col=col, bty="n")
  dev.off()
  print("Outputted to:")
  print(plot_name)
}

plot_RLE_boxes <- function(x, out_dir, species, suffix) {
  # Box plots of RLE values #
  # Normalising #
  x <- calcNormFactors(x, method = "TMM")
  x$samples$norm.factors
  
  x2 <- x
  x2$samples$norm.factors <- 1
  
  ######## Box plots of RLE values ######
  nsamples <- ncol(x)
  col <- brewer.pal(nsamples, "Paired")
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
  print("Outputted to:")
  print(plot_name)
}

plot_RLE_densities <- function(x, out_dir, species, suffix){
  # Plots RLE densities #
  
  # Normalising #
  x <- calcNormFactors(x, method = "TMM")
  x$samples$norm.factors
  
  x2 <- x
  x2$samples$norm.factors <- 1
  
  #### Density plots for RLE values ####
  nsamples <- ncol(x)
  col <- brewer.pal(nsamples, "Paired")
  prefix <- paste(out_dir, 'RLE_densities', sep='')
  plot_name <- paste(prefix, species, suffix, '.pdf', sep='_')
  pdf(plot_name)
  lcpm <- cpm(x2, log=TRUE)
  rles <- calc_cell_RLE(lcpm)
  par(mfrow=c(1,2))
  maxy=10
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
  print("Outputted to:")
  print(plot_name)
}

save_results <- function(tfit, out_dir, species, suffix, save=T) {
  # Saves the results from the tfit object #
  
  # Getting the results #
  rank_ <- order(tfit$p.value.adj)
  all_df <- as.data.frame(list('t.val'=tfit$t[rank_,1], 
                               'log2FoldChange'=tfit$coefficients[rank_,1],
                               'pvalue'=tfit$p.value[rank_,1],
                               'padj'=tfit$p.value.adj[rank_]))
  
  up_all_df <- all_df[all_df[,'t.val']>0,]
  down_all_df <- all_df[all_df[,'t.val']<0,]
  
  res_df <- all_df[all_df[,'padj']<.05,]
  if (species == 'mouse' & dim(res_df)[1]==1) { rownames(res_df) <- c(rownames(all_df)[1]) }
  up_df <- res_df[res_df[,'t.val']>0,]
  down_df <- res_df[res_df[,'t.val']<0,]
  
  # Save results for all genes #
  if (save) {
    file_name <- paste(out_dir, 'de_results_PseudoLimma', '_allGenes_', suffix, '.xlsx', sep='')
    sheet_up <- paste(species[1], '.treatment_up', sep='')
    sheet_down <- paste(species[1], '.treatment_down', sep='')
    
    write.xlsx(up_all_df, file_name, 
               sheetName=sheet_up, col.names=TRUE, row.names=TRUE, append=T)
    write.xlsx(down_all_df, file_name, 
               sheetName=sheet_down, col.names=TRUE, row.names=TRUE, append=T)
    
    # Saving results for significant genes #
    file_name <- paste(out_dir, 'de_results_PseudoLimma', '_', suffix, '.xlsx', sep='')
    sheet_up <- paste(species, '.treatment_up', sep='')
    sheet_down <- paste(species, '.treatment_down', sep='')
    
    if (nrow(up_df)!=0) {
      write.xlsx(up_df, file_name, 
               sheetName=sheet_up, col.names=TRUE, row.names=TRUE, append=T)
    }
    
    if (nrow(down_df)!=0) {
      write.xlsx(down_df, file_name, 
               sheetName=sheet_down, col.names=TRUE, row.names=TRUE, append=T)
    }
  } else {
    print(dim(up_all_df))
    print(dim(down_all_df))
  }
}





