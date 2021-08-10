# Title     : TODO
# Objective : TODO
# Created by: uqbbalde
# Created on: 21/7/21

singleR_map <- function(datas, trained, suffix,
                        out_plots, out_dir, prefix='Vlad'){
  # Performs singleR reference labelling #

  for (i in 1:length(samples)) {
    out <- classifySingleR(datas[[i]], trained, fine.tune=F, prune=T)
    
    # Checking diagnostics #
    ph <- plotScoreHeatmap(out)
    ggsave(paste0(out_plots, samples[i],'_', prefix, '_ScoreHeatmap',
                  suffix, '.pdf'), ph)
    
    # Visualising the results spatially #
    spot_names <- str_c(rep(seur_names[i], ncol(datas[[i]])), colnames(datas[[i]]))
    scores <- as.data.frame(out$scores)
    rownames(scores) <- spot_names
    labels <- out$pruned.labels
    names(labels) <- spot_names
    
    scores[,'labels'] <- labels
    write.table(scores, paste0(out_dir, samples[i], '_', prefix,
                               '_singleR_scores', suffix, '.txt'),
                quote=F, sep='\t')
  }

}




