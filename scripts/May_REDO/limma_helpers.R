# Purpose of this script is to impliment helpers for performing the Limma DE #

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


