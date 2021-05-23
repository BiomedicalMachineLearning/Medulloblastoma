# The purpose of this script is to check the FPR based on stably expressed
# genes determined by:
# https://sydneybiox.github.io/scMerge/index.html


#### Setting up environment #####
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scMerge)
  library(scater)
  library("readxl")
  library(VennDiagram)
})
source('scripts/helpers.R')

#out_dir <- 'data/DE_out/Pseudo_TMM_Limma_Voom/'
out_dir <- 'data/supps/'

#file_name <- 'data/DE_out/Pseudo_TMM_Limma_Voom/de_results_PseudoLimma.xlsx'
file_name <- 'data/supps/de_results_MDB.xlsx'

# Stable genes #
data(segList, package='scMerge')
human_stable <- segList$human$human_scSEG
mouse_stable <- segList$mouse$mouse_scSEG
lengths <- c(length(human_stable), length(mouse_stable))
max_len <- max(lengths)
rest <- rep("", max_len-min(lengths))
mouse_stable <- c(mouse_stable, rest)
stables <- as.data.frame(human_stable)
colnames(stables) <- c('human')
stables[,'mouse'] <- mouse_stable

write.table(stables,
            'data/third_party_data/stable_genes.txt', sep='\t', quote=F)

# Loading in the DE genes #
de_up <- as.data.frame(read_excel(file_name,
                                       sheet='human.treatment_up'))
rownames(de_up) <- de_up[,1]
de_up <- de_up[,-1]

de_down <- as.data.frame(read_excel(file_name,
                                       sheet='human.treatment_down'))
rownames(de_down) <- de_down[,1]
de_down <- de_down[,-1]

human_up <- getSpeciesGenes(de_up, 'hg38-', strip_prefix = T)
human_down <- getSpeciesGenes(de_down, 'hg38-', strip_prefix = T)
human_genes <- c(human_up, human_down)

mouse_up <- getSpeciesGenes(de_up, 'mm10-', strip_prefix = T)
mouse_down <- getSpeciesGenes(de_down, 'mm10-', strip_prefix = T)
mouse_genes <- c(mouse_up, mouse_down)

####### Plotting Venn diagram with stably expressed genes ########

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(segList$human$human_scSEG, segList$human$bulkRNAseqHK, human_genes),
  category.names = c("Stable" , 
                     "House-Keep" , 
                     "DE"),
  filename = paste(out_dir, 'DE_overlap_HKGenes.png', sep=''),
  output=T,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 1000, 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)




