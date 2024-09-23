library(dplyr)
library(SingleCellExperiment)
library(scran)

samples <- c("A-01_post-treatment", "O-01_pre-treatment1", "O-01_pre-treatment2", "O-01_pre-treatment3", "O-01_post-treatment", "A-03_pre-treatment", "A-03_post-treatment", "A-05_pre-treatment", "A-05_post-treatment")

k_clust <- c(7,8,9,10,11,12,13)

for(i in samples){
  
  node_meta <- read.csv(paste0("/graphsage/tmp/subgraph3NB5000/", i, "_node_meta.csv"))
  
  for(k in k_clust){
    
    trans_niche <- read.delim(paste0("/graphsage/tmp/subgraph3NB5000/gmm", k, "_5k_trained/", i, "_noMeans.txt"), header = FALSE, col.names = paste0("niche_", k))
    
    if(nrow(trans_niche) == nrow(node_meta)) {
      node_meta[, paste0("niche_", k)] <- trans_niche[, 1]
    }else{
      warning(paste("Mismatched row count in sample", i, "for k =", k))
    }
  }
  
  write.csv(node_meta, file = paste0("/analysis/output/graphsage/", i, "_niche_mapped.csv"))
  
}