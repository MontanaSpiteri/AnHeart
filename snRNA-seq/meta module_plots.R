library(msigdbr)
library(ggplot2)
library(gridExtra)
library(readxl)
library(HDF5Array)
library(scrabble)

# All code adapted from {https://github.com/parveendabas/GBMatlas/blob/master/R/Figure2.R}
# Reference: Abdelfattah et al. (2022) 'Single-cell analysis of human glioma and immune cells identifies S100A4 as an immunotherapy target' {https://www.nature.com/articles/s41467-022-28372-y}

ClusterColors <- c("#F28482", "#9D4EDD", "#6A994E", "#3A7CA5", "#EF233C",
                   "#F9C8C7", "#CDDEB9", "#A7C957", "#63A2C8", "#000000", "#FFEA00", "#F4F3EE",
                   "#FAA307", "#FFBA08", "#FFD60A", "#E85D04", "#F48C06", "#F2CC8F",
                   "#656D4A", "#A4AC86", "#B6AD90", "#7F4F24", "#A68A64", "#D6D1C1",
                   "#000000", "#ADB5BD",
                   "#9D0208")
names(ClusterColors) <- c("AC-like", "MES-like", "OPC-like", "NPC-like", "Progenitor",
                          "Astrocyte", "Oligodendrocyte", "OPC", "Neuron", "Vasculature", "Immune", "Unknown",
                          "T Cell", "NK Cell", "B Cell", "CD8 T Cell", "CD4 T Cell", "T Reg Cell", 
                          "Dendritic Cell", "Neutrophil", "Monocyte", "Microglia", "Macrophage", "Mast Cell",
                          "Endothelial", "Mural",
                          "Death")

# AnHeart Integrated dataset
sce <- loadHDF5SummarizedExperiment("/data/process/", prefix="sce_integrated_")
# remove participant A-04 (low tumor purity)
sce <- sce[, !sce$Patient == "A-04"]
sce$Patient <- factor(sce$Patient)

#Read the Genelist file #from Neftel et al. 2019
meta_modules <- read_excel("/public_datasets/neftel/metamodules.xlsx", skip = 4)
subtypelist <- lapply(meta_modules, function(column) {
  genes <- na.omit(column)
  names(genes) <- NULL
  return(genes)
})
subtypelist <- subtypelist [-c(7:8)] # filter cellcycle
subtypelist$MES <- c(subtypelist$MES1, subtypelist$MES2)
subtypelist$NPC <- c(subtypelist$NPC1, subtypelist$NPC2)
subtypelist <- subtypelist[-grep("\\d", names(subtypelist))]
rm(meta_modules)


# Across all patients - Progenitor tumor cell population only
sce <- sce[, sce$annotation %in% c("Progenitor")]
sce$annotation <- factor(sce$annotation)
sce <- sce[, !duplicated(colnames(sce))]

for(Sample_type in unique(sce$Sample_Type)){
  
  sce_sub <- sce[, sce$Sample_Type == Sample_type]
  Counts <- as.matrix(counts(sce_sub))
  CellInfo=colData(sce_sub)
  
  #get scores & filter genes also present in sce
  subtypelist <- lapply(subtypelist, function(x) {
    x[x %in% rownames(Counts)]
  })
  
  sc=score(Counts,
           groups=subtypelist,
           binmat = NULL,
           bins = NULL,
           controls = NULL,
           bin.control = F,
           center = T,
           nbin = 30,
           n = 100,
           replace = T)
  
  # save metamodule score assignments
  sc_df <- as.data.frame(sc)
  meta_module_assignments <- apply(sc_df, 1, function(row) {
    meta_module <- names(row)[which.max(row)] # Find the meta-module with the highest score
    return(meta_module)
  })
  results_df <- data.frame(
    cell_id = colnames(sce_sub),
    meta_module = meta_module_assignments
  )
  output_filename <- paste0("output/metamodule/metamodule_assignments_progenitor_", Sample_type, ".csv")
  write.csv(results_df, file = output_filename, row.names = FALSE)
  
  #get coordinates
  h=hierarchy (sc, quadrants = c("AC","OPC","MES","NPC"), log.scale = T)

  # make plots
  ##per annotation
  for( i in 1: length(names(ClusterColors))){
    name=names(ClusterColors)[i]
    color=ClusterColors[i]
    CellInfo$Clustercolor [CellInfo$annotation== name] <- color
  }
  Clustercolor <- CellInfo$Clustercolor
  names(Clustercolor ) <- row.names(CellInfo)
  colData(sce_sub)$Clustercolor <- Clustercolor
  
  #xmax=max(h$X)+(max(h$X)*0.2)
  #xmin=min(h$X)+(min(h$X)*0.2)
  #ymax=max(h$Y)+(max(h$Y)*0.6)
  #ymin=min(h$Y)+(min(h$Y)*0.6)
  
  # same axis scale for plots
  xmax <- 6
  xmin <- -3
  ymax <- 3
  ymin <- -6.5
  
  #Make the big plot--> all Annos
  groups=colData(sce_sub)[,c("Clustercolor","annotation")]
  matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
  matrix$Clustercolor <- as.character(matrix$Clustercolor)
  matrix$annotation <- as.character(matrix$annotation)
  matrix$Clustercolor[is.na(matrix$Clustercolor)] <- "gray"
  matrix$annotation[is.na(matrix$annotation)] <- "Other"
  row.names(matrix)=matrix$Row.names
  x=matrix$Clustercolor
  y=matrix$annotation
  col=x[!duplicated(x)]
  names(col)=y[!duplicated(y)]
  matrix=matrix[,-1]
  title=paste0(Sample_type)
  matrix <- as.data.frame(matrix)
  
  p0 <- ggplot(matrix, aes(x = X, y = Y, color = factor(annotation))) +
    geom_point() +
    geom_point(data = subset(matrix, annotation != "Other")) + 
    scale_color_manual(values = col, aesthetics = "colour", breaks = waiver()) +
    labs(x = "Relative meta-module score [log2(|SC1-SC2|+1)]", y = "Relative meta-module score [log2(|SC1-SC2|+1)]", colour = "annotation") +
    theme(legend.position = "none", 
          legend.text = element_text(size = 15), 
          legend.title = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(hjust = 0, vjust = -2, colour = "black", size = 10, face = "bold"),
          axis.title.y = element_text(hjust = 0, vjust = 4, colour = "black", size = 10, face = "bold"),
          panel.background = element_rect(fill = "white", colour = "white"),
          axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 2),
          plot.title = element_text(size = 25, face = "bold")) +
    scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax)) +
    scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax)) +
    ggtitle(title) +
    geom_hline(yintercept = 0, color = "black", size = 0.5) +
    geom_vline(xintercept = 0, color = "black", size = 0.5) +
    annotate("rect", xmin = xmin, xmax = xmin + 2, ymin = ymax - 1, ymax = ymax, fill = "#9D4EDD") +
    annotate("text", x = xmin + 1, y = ymax - 0.5, label = "MES-Like", color = "white", fontface = "bold", size = 4) +
    annotate("rect", xmin = xmax - 2, xmax = xmax, ymin = ymax - 1, ymax = ymax, fill = "#3A7CA5") +
    annotate("text", x = xmax - 1, y = ymax - 0.5, label = "NPC-Like", color = "white", fontface = "bold", size = 4) +
    annotate("rect", xmin = xmin, xmax = xmin + 2, ymin = ymin, ymax = ymin + 1, fill = "#F28482") +
    annotate("text", x = xmin + 1, y = ymin + 0.5, label = "AC-Like", color = "white", fontface = "bold", size = 4) +
    annotate("rect", xmin = xmax - 2, xmax = xmax, ymin = ymin, ymax = ymin + 1, fill = "#6A994E") +
    annotate("text", x = xmax - 1, y = ymin + 0.5, label = "OPC-Like", color = "white", fontface = "bold", size = 4) +
    geom_density_2d(linewidth = 0.25, colour = "black")
  P0=ggplotGrob(p0)
  
  ##make the small plots--> one per annotation
  Final <- list()
  Final[[1]]=ggplotGrob(p0)
  for (i in 1:length(unique(sce_sub$annotation))) {
    ClusterMD=colData(sce_sub)[sce_sub$annotation==paste0(unique(sce_sub$annotation)[i]),]
    groups=ClusterMD[,c("Clustercolor","annotation")]
    title=paste0(unique(sce_sub$annotation)[i])
    matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
    matrix$Clustercolor[is.na(matrix$Clustercolor)] <- "gray"
    matrix$annotation=as.character(matrix$annotation)
    matrix$annotation[is.na(matrix$annotation)] <- "Other"
    row.names(matrix)=matrix$Row.names
    x=matrix$Clustercolor
    y=matrix$annotation
    col=x[!duplicated(x)]
    names(col)=y[!duplicated(y)]
    matrix=matrix[,-1]
    matrix <- as.data.frame(matrix)
    P=ggplot(matrix, aes(x = X,
                         y =Y,color=factor(annotation)))+geom_point()+geom_point(data = subset(matrix, annotation !="Other"))+ 
      scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) + labs(x=NULL, y=NULL) + theme(legend.position = "none")+
      scale_x_continuous(expand = c(0, 0), limits = c(xmin,xmax)) + scale_y_continuous(expand = c(0, 0), limits = c(ymin,ymax))+
      theme(panel.background = element_rect(fill = "white",colour = "white"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
            axis.ticks.y=element_blank(), axis.text.y=element_blank())+
      ggtitle(title)+theme(plot.title =element_text(size=25,face="bold") )+ 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
      geom_hline(yintercept=0, color = "black", size=0.5)+
      geom_vline(xintercept=0, color = "black", size=0.5)+
      annotate("rect", xmin = xmin, xmax = xmin + 2, ymin = ymax - 1, ymax = ymax, fill = "#9D4EDD") +
      annotate("text", x = xmin + 1, y = ymax - 0.5, label = "MES-Like", color = "white", fontface = "bold", size = 1) +
      annotate("rect", xmin = xmax - 2, xmax = xmax, ymin = ymax - 1, ymax = ymax, fill = "#3A7CA5") +
      annotate("text", x = xmax - 1, y = ymax - 0.5, label = "NPC-Like", color = "white", fontface = "bold", size = 1) +
      annotate("rect", xmin = xmin, xmax = xmin + 2, ymin = ymin, ymax = ymin + 1, fill = "#F28482") +
      annotate("text", x = xmin + 1, y = ymin + 0.5, label = "AC-Like", color = "white", fontface = "bold", size = 1) +
      annotate("rect", xmin = xmax - 2, xmax = xmax, ymin = ymin, ymax = ymin + 1, fill = "#6A994E") +
      annotate("text", x = xmax - 1, y = ymin + 0.5, label = "OPC-Like", color = "white", fontface = "bold", size = 1)
    Final[[i+1]] = ggplotGrob(P)
  }
  Final[[1]]=ggplotGrob(p0)
  numofplots= length(Final)
  
  OutputDirectory <- "figures/"
  svg(file =paste0(OutputDirectory,"subtype_butterfly_progenitor_", Sample_type, ".svg"), height = 8, width =12,onefile = T)
  grid.arrange(grobs=Final, widths = c(2,1),layout_matrix = rbind(c(1, 2),c(1,3),c(1,4))) 
  dev.off()
}
