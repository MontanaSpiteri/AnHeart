library(SingleR)
library(scran)
library(scuttle)
library(SingleCellExperiment)
library(HDF5Array)

sce <- loadHDF5SummarizedExperiment("/data/process/", prefix="sce_integrated_")
sc_gbm_big <- readRDS("/public_datasets/sc_gbm_big/sce_downsampled.rds")
couturier <- readRDS("/public_datasets/Couturier_dataset_tumor_lognormalised.rds")
reference <- readRDS("/public_datasets/lister/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS")
reference$concat <- paste(reference$major_clust, reference$stage_ids, sep = "_")

# gene inputs
TopHVGs <- function(x){
  dec.parse <- modelGeneVar(x)
  TopHVGs <- getTopHVGs(dec.parse, n=5000)
  return(TopHVGs)
}
top.HVGs <- TopHVGs(sce)

ReferenceHVGs <- function(x, symbol="feature_name"){
  rownames(x) <- rowData(x)[,symbol]
  dec <- modelGeneVar(x)
  Reference_HVGs <- getTopHVGs(dec, n=5000)
  return(Reference_HVGs)
}
reference.HVGs <- ReferenceHVGs(sc_gbm_big)
reference.HVGsC <- ReferenceHVGs(couturier, symbol = "Symbol")
reference.HVGsL <- ReferenceHVGs(reference, symbol = "index")

CommonHVGs <- function(TopHVGs, ReferenceHVGs){
  CommonHVGs <- intersect(TopHVGs, ReferenceHVGs)
  return(CommonHVGs)
}
common.HVGs <- CommonHVGs(TopHVGs = top.HVGs, ReferenceHVGs = reference.HVGs)
length(common.HVGs)
common.HVGsC <- CommonHVGs(TopHVGs = top.HVGs, ReferenceHVGs = reference.HVGsC)
length(common.HVGsC)
common.HVGsL <- CommonHVGs(TopHVGs = top.HVGs, ReferenceHVGs = reference.HVGsL)
length(common.HVGsL)

# function
annotate_clust <- function(sce, reference, genes, symbol="feature_name", annotation="annotation_level_3"){
  
  #Check for excessive amount of Ensembl gene names 
  is.ens <- grepl("^ENSG", genes, ignore.case = TRUE)
  num.true <- length(is.ens[is.ens == TRUE])
  if(num.true == length(genes)){
    print(warning("Ensembl genes detected in gene names"))
  }
  
  #Check Marker Genes are actually present in dataset of interest
  is.present <- genes %in% rownames(sce)
  num.false <- length(is.present[is.present == FALSE])
  if(num.false > 0){
    print(warning("Marker genes not present in dataset of interest. Possibly wrongly named"))
  }
  
  #Annotate clusters using reference dataset
  counts(reference) <- logcounts(reference)
  rownames(reference) <- rowData(reference)[,symbol]
  available <- intersect(rownames(reference), rownames(sce))
  pred_labels <- SingleR(test=sce[available,], ref=reference[available,], restrict = genes, 
                         labels=colData(reference)[,annotation], de.method="wilcox")
  return(pred_labels$pruned.labels)
}

# run
anno.comHVGsC <- annotate_clust(sce = sce, reference = couturier, genes = common.HVGsC, symbol = "Symbol", annotation = "cluster")
save(anno.comHVGsC, file = "/data/process/anno_results_cout.RData")

anno.comHVGs <- annotate_clust(sce = sce, reference = sc_gbm_big, genes = common.HVGs)
save(anno.comHVGs, file = "/data/process/anno_results_ruiz.RData")

anno.comHVGsL <- annotate_clust(sce = sce, reference = reference, genes = common.HVGsL, symbol = "index", annotation = "concat")
save(anno.comHVGsL, file = "/data/process/anno_results_lister.RData")