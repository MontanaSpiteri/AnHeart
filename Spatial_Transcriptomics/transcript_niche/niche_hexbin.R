library(dplyr)
library(SingleCellExperiment)
library(scran)
library(ggplot2)
library(hexbin)
library(pheatmap)
library(RColorBrewer)

samples <- c("A-01_post-treatment", "O-01_pre-treatment1", "O-01_pre-treatment2", "O-01_pre-treatment3", "O-01_post-treatment", "A-03_pre-treatment", "A-03_post-treatment", "A-05_pre-treatment", "A-05_post-treatment")

#k_clust <- c(7,8,9,10,11,12,13)
# 8 clusters taken forward for subsequent analysis
k_clust <- c(8)

transcript_niche_color_list <- c(`0` = "#9FCF87",
                                 `1` = "#94FFB5", 
                                 `2` = "#9068AC",
                                 `3` = "#F77800",
                                 `4` = "#FFCC99",
                                 `5` = "#2BCE48",
                                 `6` = "#C39C94",
                                 `7` = "#005C31",
                                 `8` = "#0075DC", 
                                 `9` = "#F0A0FF", 
                                 `10` = "#C20088", 
                                 `11` = "#4C005C",
                                 `12` = "pink")


# combinations of samples and clusters
combinations <- expand.grid(sample = samples, k = k_clust, stringsAsFactors = FALSE)


plot_sample_hexbin <- function(sample, k, binwidth=70) {
  
  input_dir <- "/analysis/output/graphsage/"
  dir <- paste0(input_dir, "niche_", k, "_binwidth", binwidth, "/")
  
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  
  node_meta <- read.csv(paste0(input_dir, sample, "_niche_mapped.csv"))
  
  # number of bins required to achieve consistent bin size relative to data range
  x_range <- range(node_meta$x_location)
  xbins <- round((x_range[2] - x_range[1]) / binwidth)
  
  bins <- hexbin(x = node_meta$x_location, y = node_meta$y_location, xbins = xbins, IDs = TRUE)
  node_meta$bin <- bins@cID
  
  # majority label for each bin
  majority_label <- function(labels) {
    names(sort(table(labels), decreasing = TRUE))[1]
  }
  bin_labels <- tapply(node_meta[paste0("niche_", k)], node_meta$bin, majority_label)
  bin_labels <- data.frame(bin = names(bin_labels), label = bin_labels, stringsAsFactors = FALSE)
  
  bin_centers <- data.frame(hcell2xy(bins), count = bins@count, bin = bins@cell)
  plot_data <- merge(bin_centers, bin_labels, by = "bin")
  plot_data$niche <- ifelse(plot_data$count < 10, "white", as.character(plot_data$label)) # filter bins with less than 10 transcripts
  
  # Calculate and save niche percentages
  niche_percentage <- plot_data %>%
    group_by(niche) %>%
    summarise(count = n()) %>% 
    ungroup() %>%
    mutate(total_bins = sum(count),
           percentage = (count / total_bins) * 100)
  write.csv(niche_percentage, file = paste0(dir, sample, "_niche_percentage.csv"), row.names = FALSE)
  
  # save hexbin info
  write.csv(plot_data, file = paste0(dir, sample, "_hexbin_info.csv"), col.names = TRUE)
  
  plot <- ggplot(plot_data, aes(x = x, y = y, fill = factor(niche))) +
    geom_hex(stat = "identity") +
    scale_fill_manual(values = c("white" = "white", transcript_niche_color_list)) +
    theme_classic() +
    labs(fill = "Transcript Niche") + 
    theme(legend.position = "right")
  
  # save hexbin transcript assignments
  write.csv(node_meta$bin, file = paste0(dir, sample, "_hexbins.csv"))
  
  ggsave(paste0(dir, sample, "_hexbin_plot.png"), plot = plot)
  
  gc()
  
}

mapply(plot_sample_hexbin, combinations$sample, combinations$k)